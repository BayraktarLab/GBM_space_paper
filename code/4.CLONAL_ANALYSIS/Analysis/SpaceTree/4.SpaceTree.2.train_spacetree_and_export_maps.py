"""Train SpaceTree models and export Visium clone/cell-type maps.

Inputs
------
data/processed/data_<sample>.pt
data/processed/full_encoding_<sample>.pkl
data/raw/visium_data/*/*/filtered_feature_bc_matrix.h5

Outputs
-------
notebooks/results/<section>_<sample>_clones.png
notebooks/results/<section>_<sample>_cell_types.png
data/processed/visium_annotated_<sample>.h5ad
"""

from __future__ import annotations

import argparse
import glob
import pickle
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import torch
import lightning.pytorch as pl
from lightning.pytorch.loggers import TensorBoardLogger
from sklearn.model_selection import train_test_split
from torch_geometric.loader import NeighborLoader

import spaceTree.utils as utils
from pipeline_config import PROCESSED_DIR, RAW_DIR, RESULTS_DIR, SKIP_VISIUM_SAMPLES, parse_patients, section_short
from spaceTree.models import GATLightningModule_sampler


warnings.simplefilter("ignore")


def compute_class_weights(y_train):
    y_train = y_train.detach().cpu().numpy() if torch.is_tensor(y_train) else np.asarray(y_train)
    class_sample_count = np.array([len(np.where(y_train == value)[0]) for value in np.unique(y_train)])
    return 1.0 / class_sample_count


def balanced_split(data, hold_in, size=0.3):
    y_type = data.y_type.detach().cpu().numpy()
    y_clone = data.y_clone.detach().cpu().numpy()
    train_indices_type, test_indices_type, _, _ = train_test_split(
        hold_in,
        y_type[hold_in],
        test_size=0.5,
        stratify=data.y_type[hold_in],
        random_state=42,
    )

    train_indices_final, test_indices_final = [], []
    for subset_indices in [train_indices_type, test_indices_type]:
        train_subset, test_subset, _, _ = train_test_split(
            subset_indices,
            y_clone[subset_indices],
            test_size=size,
            stratify=data.y_clone[subset_indices],
            random_state=42,
        )
        train_indices_final.extend(train_subset)
        test_indices_final.extend(test_subset)
    return train_indices_final, test_indices_final


def check_class_distributions(data, weight_clone, weight_type, norm_sim, no_diploid=False):
    num_clone_train = data.y_clone[data.train_mask].unique().shape[0]
    num_clone_total = len(data.y_clone.unique())
    assert num_clone_total - 1 == num_clone_train, (
        f"Clone classes in training set ({num_clone_train}) do not match total labelled "
        f"clone classes ({num_clone_total - 1})"
    )
    assert num_clone_total - 1 == len(weight_clone), "Clone weights do not match clone classes"
    if no_diploid:
        assert num_clone_total - 1 == norm_sim.shape[0] - 1, "Similarity matrix shape mismatch"
    else:
        assert num_clone_total - 1 == norm_sim.shape[0], "Similarity matrix shape mismatch"

    num_type_train = data.y_type[data.train_mask].unique().shape[0]
    num_type_total = len(data.y_type.unique())
    assert num_type_total - 1 == num_type_train, "Cell-type classes in training set mismatch"


def parse_visium_path(path: str) -> tuple[str, str]:
    parts = path.rstrip("/").split("/")
    return parts[-2], parts[-1]


def load_visium(path: str, sample: str) -> sc.AnnData:
    visium = sc.read_visium(
        path,
        genome=None,
        count_file="filtered_feature_bc_matrix.h5",
        library_id=None,
        load_images=True,
        source_image_path=None,
    )
    visium.var_names_make_unique()
    visium.obsm["spatial"] = np.asarray(visium.obsm["spatial"], dtype=int)
    visium.obs.index = [f"{sample}_{barcode}" for barcode in visium.obs.index]
    return visium


def choose_device(gpu: int | None) -> tuple[torch.device, str, list[int] | str]:
    if torch.cuda.is_available() and gpu is not None:
        return torch.device(f"cuda:{gpu}"), "gpu", [gpu]
    if torch.cuda.is_available():
        return torch.device("cuda:0"), "gpu", [0]
    return torch.device("cpu"), "cpu", "auto"


def train_and_predict(path: str, args) -> None:
    section, sample = parse_visium_path(path)
    if sample in SKIP_VISIUM_SAMPLES:
        print(f"[train] {section}: skipping excluded sample {sample}")
        return

    clones_png = RESULTS_DIR / f"{section}_{sample}_clones.png"
    cell_types_png = RESULTS_DIR / f"{section}_{sample}_cell_types.png"
    annotated_out = PROCESSED_DIR / f"visium_annotated_{sample}.h5ad"
    if clones_png.exists() and cell_types_png.exists() and annotated_out.exists() and not args.force:
        print(f"[train] {section}: outputs exist for {sample}; use --force to rebuild")
        return

    data_path = PROCESSED_DIR / f"data_{sample}.pt"
    encoding_path = PROCESSED_DIR / f"full_encoding_{sample}.pkl"
    if not data_path.exists() or not encoding_path.exists():
        print(f"[train] {section}: missing graph inputs for {sample}")
        return

    print(f"[train] {section}: training model for {sample}")
    visium_full = load_visium(path, sample)
    try:
        data = torch.load(data_path, weights_only=False)
    except TypeError:
        data = torch.load(data_path)
    with open(encoding_path, "rb") as handle:
        encoder_dict = pickle.load(handle)

    node_encoder_rev = {val: key for key, val in encoder_dict["nodes"].items()}
    node_encoder_clone = {val: key for key, val in encoder_dict["clones"].items()}
    node_encoder_ct = {val: key for key, val in encoder_dict["types"].items()}

    data.edge_attr = data.edge_attr.reshape((-1, 1))
    hold_out_indices = np.where(data.y_clone == -1)[0]
    hold_out = torch.tensor(hold_out_indices, dtype=torch.long)
    hold_in = [idx for idx in np.arange(data.x.shape[0]) if idx not in hold_out_indices]

    train_indices, test_indices = balanced_split(data, hold_in, size=args.test_size)
    data.train_mask = torch.tensor(train_indices, dtype=torch.long)
    data.test_mask = torch.tensor(test_indices, dtype=torch.long)
    data.hold_out = hold_out

    weight_type = torch.tensor(compute_class_weights(data.y_type[data.train_mask]), dtype=torch.float)
    weight_clone = torch.tensor(compute_class_weights(data.y_clone[data.train_mask]), dtype=torch.float)
    data.num_classes_clone = len(data.y_clone.unique())
    data.num_classes_type = len(data.y_type.unique())
    if hasattr(data, "edge_type"):
        del data.edge_type

    train_loader = NeighborLoader(data, num_neighbors=[10] * 3, batch_size=args.batch_size, input_nodes=data.train_mask)
    valid_loader = NeighborLoader(data, num_neighbors=[10] * 3, batch_size=args.batch_size, input_nodes=data.test_mask)

    norm_sim = torch.tensor(np.diag(np.ones(data.num_classes_clone - 1)), dtype=torch.float)
    device, accelerator, devices = choose_device(args.gpu)
    data = data.to(device)
    weight_clone = weight_clone.to(device)
    weight_type = weight_type.to(device)
    norm_sim = norm_sim.to(device)

    diploid_indices = [key for key, val in node_encoder_clone.items() if val == "diploid"]
    test_clone_labels = data.y_clone[test_indices].detach().cpu().unique().tolist()
    no_diploid = bool(diploid_indices) and diploid_indices[0] not in test_clone_labels
    check_class_distributions(data, weight_clone, weight_type, norm_sim, no_diploid=no_diploid)

    model = GATLightningModule_sampler(
        data,
        weight_clone,
        weight_type,
        norm_sim=norm_sim,
        learning_rate=args.learning_rate,
        heads=args.heads,
        dim_h=args.hidden_dim,
    ).to(device)
    logger = TensorBoardLogger("logs_visium", name=f"round1_{section}_{sample}")
    early_stop = pl.callbacks.EarlyStopping(
        monitor="validation_combined_loss",
        min_delta=1e-4,
        patience=args.patience,
        verbose=True,
        mode="min",
    )
    trainer = pl.Trainer(
        max_epochs=args.max_epochs,
        accelerator=accelerator,
        devices=devices,
        logger=logger,
        callbacks=[early_stop],
        log_every_n_steps=10,
    )
    trainer.fit(model, train_loader, valid_loader)

    model.eval()
    with torch.no_grad():
        out, _, _ = model(data)

    clone_res, ct_res = utils.get_results(
        out, data, node_encoder_rev, node_encoder_ct, node_encoder_clone, activation="softmax"
    )
    visium = visium_full[~visium_full.obs.join(clone_res)[clone_res.columns[0]].isna()].copy()
    visium.obs = visium.obs.join(clone_res).join(ct_res)
    visium.obs["clone"] = visium.obs[clone_res.columns].idxmax(axis=1).astype(str)
    visium.obs["cell_type"] = visium.obs[ct_res.columns].idxmax(axis=1).astype(str)
    visium = visium[~visium.obs.cell_type.isna()].copy()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    with mpl.rc_context({"axes.facecolor": "black", "figure.figsize": [10, 10], "font.size": 30}):
        sc.pl.spatial(visium, color=clone_res.columns, wspace=0.3, img_key="lowres", alpha_img=0.5, show=False)
        plt.savefig(clones_png, bbox_inches="tight", dpi=180)
        plt.close("all")

    with mpl.rc_context({"axes.facecolor": "black", "figure.figsize": [10, 10], "font.size": 30}):
        sc.pl.spatial(visium, color=ct_res.columns, wspace=0.3, img_key="lowres", alpha_img=0.5, show=False)
        plt.savefig(cell_types_png, bbox_inches="tight", dpi=180)
        plt.close("all")

    visium.obs.columns = [str(x) for x in visium.obs.columns]
    visium.obs["section"] = section_short(section)
    visium.write_h5ad(annotated_out)
    del model
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    print(f"[train] {section}: wrote plots and {annotated_out}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--patients", default="all", help="Comma-separated patients or 'all'.")
    parser.add_argument("--gpu", type=int, default=None, help="GPU index; defaults to cuda:0 if available.")
    parser.add_argument("--max-epochs", type=int, default=1000)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--learning-rate", type=float, default=5e-3)
    parser.add_argument("--hidden-dim", type=int, default=200)
    parser.add_argument("--heads", type=int, default=1)
    parser.add_argument("--test-size", type=float, default=0.3)
    parser.add_argument("--force", action="store_true", help="Rebuild outputs even if they exist.")
    parser.add_argument("--strict", action="store_true", help="Stop at the first failed sample.")
    args = parser.parse_args()

    visium_paths = sorted(glob.glob(str(RAW_DIR / "visium_data" / "*" / "*" / "")))
    for patient in parse_patients(args.patients):
        patient_paths = [path for path in visium_paths if parse_visium_path(path)[0].startswith(patient)]
        for path in patient_paths:
            try:
                train_and_predict(path, args)
            except Exception as exc:
                message = f"[train] failed for {path}: {exc}"
                if args.strict:
                    raise RuntimeError(message) from exc
                print(message)


if __name__ == "__main__":
    main()
