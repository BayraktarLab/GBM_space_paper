"""Match RNA and ATAC clone calls into joint clone labels.

Inputs
------
data/interim/<patient_lower>_annotation_5k.h5ad
data/interim/atac_annotated_<PATIENT>.h5ad

Outputs
-------
data/interim/<patient_lower>_joint_annotation2.h5ad
notebooks/results/joint_<PATIENT>_*.png
"""

from __future__ import annotations

import argparse
import warnings

import infercnvpy as cnv
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.interpolate import interp1d
from scipy.spatial import distance
from sklearn.semi_supervised import LabelSpreading

from pipeline_config import INTERIM_DIR, RESULTS_DIR, parse_patients


warnings.simplefilter("ignore")
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def as_array(matrix) -> np.ndarray:
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


def save_current_figure(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight", dpi=180)
    plt.close("all")


def score_mat(
    x_rna: pd.DataFrame,
    x_atac: pd.DataFrame,
    adata_rna: sc.AnnData,
    adata_atac: sc.AnnData,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    overlap_mat = np.zeros((x_rna.shape[0], x_atac.shape[0]))
    cosine_mat = np.zeros((x_rna.shape[0], x_atac.shape[0]))
    for i, cluster_i in enumerate(x_rna.index):
        for j, cluster_j in enumerate(x_atac.index):
            cells_i = adata_rna.obs.index[adata_rna.obs.rna_clone == cluster_i]
            cells_j = adata_atac.obs.index[adata_atac.obs.atac_clone == cluster_j]
            denom = max(1, min(len(cells_i), len(cells_j)))
            overlap_mat[i, j] = len(set(cells_i).intersection(cells_j)) / denom
            cosine = 1 - distance.cosine(x_rna.loc[cluster_i], x_atac.loc[cluster_j])
            cosine_mat[i, j] = 0 if np.isnan(cosine) else cosine
    return (
        pd.DataFrame(overlap_mat, index=x_rna.index, columns=x_atac.index),
        pd.DataFrame(cosine_mat, index=x_rna.index, columns=x_atac.index),
    )


def build_clone_centroids(
    adata_rna: sc.AnnData, adata_atac: sc.AnnData, shared_cells: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    x_rna_raw = as_array(adata_rna.obsm["X_cnv"])
    rows, cols = x_rna_raw.shape
    new_cols = adata_atac.X.shape[1]
    old_col_indices = np.linspace(0, cols - 1, cols)
    new_col_indices = np.linspace(0, cols - 1, new_cols)
    new_data = np.zeros((rows, new_cols))
    for i in range(rows):
        interpolate_function = interp1d(old_col_indices, x_rna_raw[i, :], kind="linear")
        new_data[i, :] = interpolate_function(new_col_indices)

    x_rna = pd.DataFrame(new_data, index=adata_rna.obs.index)
    x_rna = x_rna.loc[shared_cells].join(adata_rna.obs["rna_clone"])
    x_rna = x_rna.groupby("rna_clone").mean()

    x_atac = pd.DataFrame(as_array(adata_atac.X), index=adata_atac.obs.index)
    x_atac = x_atac.loc[shared_cells].join(adata_atac.obs["atac_clone"])
    x_atac = x_atac.groupby("atac_clone").mean()
    return x_rna, x_atac


def merge_similar_rna_clones(score: pd.DataFrame, min_distance: float) -> dict[str, str]:
    if score.shape[0] < 2:
        return {}
    rna_clone_dist = distance.squareform(distance.pdist(score, metric="correlation"))
    rna_clone_dist = pd.DataFrame(rna_clone_dist, index=score.index, columns=score.index)
    pairs = []
    for i in range(len(rna_clone_dist)):
        for j in range(i + 1, len(rna_clone_dist)):
            if rna_clone_dist.iloc[i, j] < min_distance:
                pairs.append({rna_clone_dist.index[i], rna_clone_dist.columns[j]})

    groups = []
    while pairs:
        first, *rest = pairs
        first = set(first)
        changed = True
        while changed:
            changed = False
            for other in rest[:]:
                if first.intersection(other):
                    first |= other
                    rest.remove(other)
                    changed = True
        groups.append(first)
        pairs = rest

    replace = {}
    next_id = score.shape[0]
    for group in groups:
        new_label = str(next_id)
        for clone in group:
            replace[clone] = new_label
        next_id += 1
    return replace


def match_clones(score: pd.DataFrame, quantile: float) -> dict[str, list[str]]:
    melt = score.reset_index().melt(id_vars="rna_clone", var_name="atac_clone", value_name="value")
    threshold = np.quantile(melt.value, quantile)
    melt = melt[melt.value > threshold].sort_values(by="value", ascending=False)
    matches: dict[str, list[str]] = {}
    atac_matched = set()
    for tup in melt.itertuples(index=False):
        if tup.atac_clone in atac_matched:
            continue
        matches.setdefault(tup.rna_clone, []).append(tup.atac_clone)
        atac_matched.add(tup.atac_clone)
    return matches


def assign_seed_joint_clones(
    adata_rna: sc.AnnData, adata_atac: sc.AnnData, matches: dict[str, list[str]]
) -> pd.DataFrame:
    new_clone_dict = {}
    new_id = 0
    for rna_clone, atac_clones in matches.items():
        rna_cells = set(adata_rna.obs.index[adata_rna.obs.rna_clone == rna_clone])
        atac_cells = set(adata_atac.obs.index[adata_atac.obs.atac_clone.isin(atac_clones)])
        joint_cells = rna_cells.intersection(atac_cells)
        clone_id = "diploid" if rna_clone == "diploid" else str(new_id)
        for cell in joint_cells:
            new_clone_dict[cell] = clone_id
        if rna_clone != "diploid":
            new_id += 1
    return pd.DataFrame.from_dict(new_clone_dict, orient="index", columns=["joint_clone"])


def propagate_joint_clones(adata_rna: sc.AnnData) -> sc.AnnData:
    x = as_array(adata_rna.obsm["X_cnv"])
    y = adata_rna.obs["joint_clone"].astype(str).replace("nan", "-1")
    labels = sorted([label for label in y.unique() if label != "-1"])
    label_mapping = {label: i for i, label in enumerate(labels)}
    label_mapping["-1"] = -1
    y_encoded = y.replace(label_mapping).astype(int)

    model = LabelSpreading()
    model.fit(x, y_encoded)
    reversed_labels = {v: k for k, v in label_mapping.items()}
    predicted = [reversed_labels[i] for i in model.transduction_]
    probabilities = model.label_distributions_[
        np.arange(len(model.transduction_)), model.transduction_
    ]
    adata_rna.obs["predicted"] = predicted
    adata_rna.obs["predicted_probability"] = probabilities
    return adata_rna


def plot_joint_qc(
    adata_rna: sc.AnnData,
    adata_atac: sc.AnnData,
    score: pd.DataFrame,
    patient: str,
) -> None:
    sns.clustermap(score, cmap="Blues")
    save_current_figure(RESULTS_DIR / f"joint_{patient}_rna_atac_score.png")

    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
    cnv.pl.umap(adata_rna[~adata_rna.obs.joint_clone.isna()], color=["joint_clone"], s=20, ax=axes[0], show=False)
    cnv.pl.umap(adata_rna, color=["predicted"], s=2, ax=axes[1], show=False)
    cnv.pl.umap(adata_rna, color=["predicted_probability"], s=2, ax=axes[2], show=False)
    save_current_figure(RESULTS_DIR / f"joint_{patient}_label_propagation.png")

    if adata_atac.n_obs:
        cnv.pl.chromosome_heatmap(adata_atac, groupby="atac_clone", dendrogram=True, figsize=(16, 10), show=False)
        save_current_figure(RESULTS_DIR / f"joint_{patient}_atac_clone_heatmap.png")


def process_patient(patient: str, alpha: float, quantile: float, min_distance: float) -> None:
    rna_path = INTERIM_DIR / f"{patient.lower()}_annotation_5k.h5ad"
    atac_path = INTERIM_DIR / f"atac_annotated_{patient}.h5ad"
    if not rna_path.exists():
        raise FileNotFoundError(f"Missing RNA annotation for {patient}: {rna_path}")
    if not atac_path.exists():
        raise FileNotFoundError(f"Missing ATAC annotation for {patient}: {atac_path}")

    print(f"[joint] {patient}: reading RNA and ATAC annotations")
    adata_rna = sc.read_h5ad(rna_path)
    adata_rna.obs.index = adata_rna.obs.cell_id.astype(str)
    adata_rna.obs.rename(columns={"cnv_leiden": "rna_clone"}, inplace=True)

    adata_atac = sc.read_h5ad(atac_path)
    adata_atac.obs.rename(columns={"cnv_leiden": "atac_clone"}, inplace=True)
    shared_cells = sorted(set(adata_rna.obs.index).intersection(adata_atac.obs.index))
    if not shared_cells:
        raise ValueError(f"No shared RNA/ATAC cell ids for {patient}")

    adata_rna = adata_rna[shared_cells].copy()
    adata_atac = adata_atac[shared_cells].copy()
    adata_atac.uns["cnv"] = {
        "chr_pos": adata_atac.var.groupby("chromosome")["pos"].min().to_dict()
    }

    x_rna, x_atac = build_clone_centroids(adata_rna, adata_atac, shared_cells)
    overlap_mat, cosine_mat = score_mat(x_rna, x_atac, adata_rna, adata_atac)
    score = alpha * overlap_mat + (1 - alpha) * cosine_mat

    replacements = merge_similar_rna_clones(score, min_distance)
    if replacements:
        adata_rna.obs["rna_clone_original"] = adata_rna.obs["rna_clone"]
        adata_rna.obs["rna_clone"] = adata_rna.obs["rna_clone"].replace(replacements)
        x_rna, x_atac = build_clone_centroids(adata_rna, adata_atac, shared_cells)
        overlap_mat, cosine_mat = score_mat(x_rna, x_atac, adata_rna, adata_atac)
        score = alpha * overlap_mat + (1 - alpha) * cosine_mat

    matches = match_clones(score, quantile)
    if not matches:
        raise ValueError(f"No RNA/ATAC clone matches passed the score quantile for {patient}")
    new_clones = assign_seed_joint_clones(adata_rna, adata_atac, matches)
    adata_rna.obs = adata_rna.obs.join(new_clones)
    adata_rna = propagate_joint_clones(adata_rna)
    adata_atac.obs = adata_atac.obs.join(new_clones)
    plot_joint_qc(adata_rna, adata_atac, score, patient)

    adata_rna.obs["joint_clone"] = adata_rna.obs["joint_clone"].astype(object).fillna("unmatched")
    for column in adata_rna.obs.select_dtypes(include=["category"]).columns:
        adata_rna.obs[column] = adata_rna.obs[column].astype(str)

    out = INTERIM_DIR / f"{patient.lower()}_joint_annotation2.h5ad"
    adata_rna.write_h5ad(out)
    print(f"[joint] {patient}: wrote {out}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--patients", default="all", help="Comma-separated patients or 'all'.")
    parser.add_argument("--alpha", type=float, default=0.5, help="Weight of cell-overlap score.")
    parser.add_argument("--match-quantile", type=float, default=0.8)
    parser.add_argument("--merge-distance", type=float, default=0.01)
    args = parser.parse_args()

    for patient in parse_patients(args.patients):
        process_patient(patient, args.alpha, args.match_quantile, args.merge_distance)


if __name__ == "__main__":
    main()
