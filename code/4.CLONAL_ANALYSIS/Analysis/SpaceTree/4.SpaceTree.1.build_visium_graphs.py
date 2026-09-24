"""Build joint scRNA-to-Visium graph datasets for all samples.

Inputs
------
data/interim/<patient_lower>_joint_annotation2.h5ad
data/raw/visium_data/*/*/filtered_feature_bc_matrix.h5
data/interim/res_scvi_<sample>_new.csv, or enough data to create it with
spaceTree.preprocessing.run_scvi

Outputs
-------
data/interim/embedding_spatial_<sample>.csv
data/interim/embedding_rna_<sample>.csv
data/processed/data_<sample>.pt
data/processed/full_encoding_<sample>.pkl
"""

from __future__ import annotations

import argparse
import glob
import pickle
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from torch_geometric.loader import NeighborLoader

import spaceTree.dataset as dataset
import spaceTree.preprocessing as pp
from pipeline_config import (
    INTERIM_DIR,
    PROCESSED_DIR,
    RAW_DIR,
    SECTION_RENAMES,
    SKIP_VISIUM_SAMPLES,
    parse_patients,
    section_short,
)


warnings.simplefilter("ignore")


def parse_visium_path(path: str) -> tuple[str, str]:
    parts = path.rstrip("/").split("/")
    return parts[-2], parts[-1]


def load_joint_annotation(patient: str, min_probability: float) -> sc.AnnData:
    path = INTERIM_DIR / f"{patient.lower()}_joint_annotation2.h5ad"
    if not path.exists():
        raise FileNotFoundError(f"Missing joint RNA annotation for {patient}: {path}")

    adata = sc.read_h5ad(path)
    adata.obs["source"] = "scRNA"
    if "TME_coarse_GBM_granular" in adata.obs:
        adata.obs["cell_type"] = adata.obs["TME_coarse_GBM_granular"]
    elif "cell_type" not in adata.obs:
        raise ValueError(f"{path} has no cell type column")

    if "predicted" in adata.obs:
        adata.obs["cnv_leiden"] = adata.obs["predicted"]
    elif "cnv_leiden" not in adata.obs:
        raise ValueError(f"{path} has no predicted/cnv_leiden clone column")

    if "predicted_probability" in adata.obs:
        adata = adata[adata.obs.predicted_probability.astype(float) >= min_probability].copy()

    adata.obs["section"] = adata.obs["section"].astype(str).replace(SECTION_RENAMES)
    return adata


def load_visium(path: str, section: str, sample: str) -> sc.AnnData:
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
    visium.obs["source"] = "visium"
    visium.var["mt"] = visium.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(visium, qc_vars=["mt"], inplace=True)
    sc.pp.filter_cells(visium, min_counts=2000)
    sc.pp.filter_genes(visium, min_cells=10)
    visium.obs.index = [f"{sample}_{barcode}" for barcode in visium.obs.index]
    visium.obs["section"] = section
    visium.obs["sample"] = sample
    return visium


def run_or_load_scvi(combined: sc.AnnData, sample: str) -> pd.DataFrame:
    path = INTERIM_DIR / f"res_scvi_{sample}_new.csv"
    if path.exists():
        return pd.read_csv(path, index_col=0)
    return pp.run_scvi(combined, str(path))


def build_dataset_for_visium(
    adata_seq: sc.AnnData,
    visium_path: str,
    force: bool,
) -> None:
    section, sample = parse_visium_path(visium_path)
    if sample in SKIP_VISIUM_SAMPLES:
        print(f"[prep] {section}: skipping excluded sample {sample}")
        return

    data_out = PROCESSED_DIR / f"data_{sample}.pt"
    encoding_out = PROCESSED_DIR / f"full_encoding_{sample}.pkl"
    if data_out.exists() and encoding_out.exists() and not force:
        print(f"[prep] {section}: outputs exist for {sample}; use --force to rebuild")
        return

    short = section_short(section)
    if short not in set(adata_seq.obs.section):
        print(f"[prep] {section}: no matching scRNA section {short}")
        return

    print(f"[prep] {section}: building graph for sample {sample}")
    sc_section = adata_seq[adata_seq.obs.section == short].copy()
    sc_section.X = sc_section.layers["counts"]
    sc_section.var["col"] = "c"
    sc_section.obs = sc_section.obs[["source", "section", "cnv_leiden", "cell_type"]].copy()

    visium = load_visium(visium_path, section, sample)
    combined = sc_section.concatenate(visium)
    cell_source = run_or_load_scvi(combined, sample)
    cell_source.index = [x[:-2] for x in cell_source.index]

    emb_spatial = cell_source[cell_source.source == "visium"][["0", "1"]]
    emb_spatial.index = [f"{sample}_{x.split('_')[1]}" for x in emb_spatial.index]
    emb_rna = cell_source[cell_source.source == "scRNA"][["0", "1"]]
    emb_rna.index = [f"{x}_{sample}" for x in emb_rna.index]

    edges_sc2vis = pp.record_edges(emb_spatial, emb_rna, 10, "sc2vis")
    try:
        pp.show_weights_distribution(edges_sc2vis, visium, "visium")
    except ValueError:
        visium.obs = visium.obs.drop(columns=["weight"], errors="ignore")
        pp.show_weights_distribution(edges_sc2vis, visium, "visium")
    edges_sc2vis = pp.normalize_edge_weights(edges_sc2vis)

    edges_sc2sc = pp.normalize_edge_weights(pp.record_edges(emb_rna, emb_rna, 10, "sc2sc"))
    visium.obs.array_row = visium.obs.array_row.astype(int)
    visium.obs.array_col = visium.obs.array_col.astype(int)
    edges_vis2grid = pp.create_edges_for_visium_nodes(visium)
    edges = pd.concat([edges_sc2vis, edges_sc2sc, edges_vis2grid])
    edges.node1 = edges.node1.astype(str)
    edges.node2 = edges.node2.astype(str)

    pp.save_edges_and_embeddings(edges, emb_spatial, emb_rna, suffix=sample, outdir=str(INTERIM_DIR))

    overcl_source = sc_section.obs[["cnv_leiden", "cell_type", "section"]].copy()
    overcl_source = overcl_source.reset_index()
    overcl_source.columns = ["node1", "clone", "cell_type", "section"]

    tmp = emb_rna.copy()
    tmp["node1"] = ["_".join(x.split("_")[:3]) for x in tmp.index]
    tmp = tmp.reset_index().merge(overcl_source, on="node1", how="left")[
        ["index", "clone", "cell_type", "section"]
    ]
    tmp.columns = ["node1", "clone", "cell_type", "section"]

    edges, overcl = dataset.preprocess_data(edges, tmp, "sc2vis", "vis2grid")
    embedding_paths = {
        "spatial": str(INTERIM_DIR / f"embedding_spatial_{sample}.csv"),
        "rna": str(INTERIM_DIR / f"embedding_rna_{sample}.csv"),
    }
    edges = edges[~edges.node2.isin(set(edges.node2).difference(set(edges.node1)))]
    emb_vis_nodes, emb_rna_nodes, edges, node_encoder = dataset.read_and_merge_embeddings(
        embedding_paths, edges, 10
    )

    edges.weight = edges.weight.astype(float)
    edges.clone = edges.clone.astype(str).replace("nan", np.nan)
    edges.cell_type = edges.cell_type.astype(str).replace("nan", np.nan)

    data, encoding_dict = dataset.create_data_object(edges, emb_vis_nodes, emb_rna_nodes, node_encoder)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    torch.save(data, data_out)
    with open(encoding_out, "wb") as fp:
        pickle.dump(encoding_dict, fp)

    data.edge_attr = data.edge_attr.reshape((-1, 1))
    hold_out_indices = np.where(data.y_clone == -1)[0]
    hold_in = [idx for idx in np.arange(data.x.shape[0]) if idx not in hold_out_indices]
    if hasattr(data, "edge_type"):
        del data.edge_type
    loader = NeighborLoader(data, num_neighbors=[10] * 3, batch_size=128, input_nodes=hold_in)
    for batch in loader:
        assert -1 not in batch.y_clone.unique()
    print(f"[prep] {section}: wrote {data_out}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--patients", default="all", help="Comma-separated patients or 'all'.")
    parser.add_argument("--min-probability", type=float, default=0.9)
    parser.add_argument("--force", action="store_true", help="Rebuild outputs even if they exist.")
    parser.add_argument("--strict", action="store_true", help="Stop at the first failed sample.")
    args = parser.parse_args()

    visium_paths = sorted(glob.glob(str(RAW_DIR / "visium_data" / "*" / "*" / "")))
    if not visium_paths:
        raise FileNotFoundError("No Visium sections found under data/raw/visium_data/*/*/")

    for patient in parse_patients(args.patients):
        print(f"[prep] {patient}: loading joint RNA annotation")
        adata_seq = load_joint_annotation(patient, args.min_probability)
        patient_paths = [path for path in visium_paths if parse_visium_path(path)[0].startswith(patient)]
        for path in patient_paths:
            try:
                build_dataset_for_visium(adata_seq, path, args.force)
            except Exception as exc:
                message = f"[prep] failed for {path}: {exc}"
                if args.strict:
                    raise RuntimeError(message) from exc
                print(message)


if __name__ == "__main__":
    main()
