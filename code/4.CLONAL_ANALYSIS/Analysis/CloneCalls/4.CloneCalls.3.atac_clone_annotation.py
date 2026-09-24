"""Annotate ATAC CNV-like profiles for all GBM patients.

Inputs
------
data/interim/atac_<PATIENT>.h5ad
data/interim/annotation_obs.29_06_23.csv

Outputs
-------
data/interim/atac_annotated_<PATIENT>.h5ad
notebooks/results/atac_<PATIENT>_*.png
"""

from __future__ import annotations

import argparse
import warnings

import infercnvpy as cnv
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc

from pipeline_config import INTERIM_DIR, RESULTS_DIR, parse_patients


warnings.simplefilter("ignore")
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_annotation() -> pd.DataFrame:
    anno = pd.read_csv(INTERIM_DIR / "annotation_obs.29_06_23.csv", index_col=0)
    anno["cell_type"] = anno["TME_coarse_GBM_granular"].astype(str).str.split(".").str[0]
    return anno


def add_cnv_metadata(adata: sc.AnnData) -> sc.AnnData:
    adata.var.columns = ["chromosome", "start", "end"]
    adata.var["pos"] = np.arange(adata.var.shape[0])
    adata.obsm["X_cnv"] = adata.X
    adata.uns["cnv"] = {
        "chr_pos": adata.var.groupby("chromosome")["pos"].min().to_dict()
    }
    return adata


def save_current_figure(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight", dpi=180)
    plt.close("all")


def plot_qc(adata: sc.AnnData, patient: str) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    sc.pl.umap(adata, color="cnv_score", ax=axes[0], show=False)
    sc.pl.umap(adata, color="cell_type", ax=axes[1], show=False)
    save_current_figure(RESULTS_DIR / f"atac_{patient}_umap_cell_type.png")

    cnv.pl.chromosome_heatmap_summary(
        adata, groupby="cnv_leiden", dendrogram=True, figsize=(16, 10), show=False
    )
    save_current_figure(RESULTS_DIR / f"atac_{patient}_cnv_summary.png")

    df = adata.obs[["cell_type", "cnv_leiden", "cnv_score"]].copy()
    df.columns = ["cell_type", "clone", "cnv_score"]
    grouped = df.groupby("clone")["cell_type"].value_counts(normalize=True).unstack().fillna(0)
    cnv_score = df.groupby("clone")["cnv_score"].mean().sort_values(ascending=False)

    ax = grouped.loc[cnv_score.index].plot(kind="bar", stacked=True, figsize=(10, 7))
    ax.plot(range(len(cnv_score)), cnv_score.values, label="mean cnv_score", ls="--", color="black")
    ax.axhline(0.05, label="diploid cutoff", ls="--", color="red")
    ax.set_ylabel("Fraction of cells / CNV score")
    ax.legend(bbox_to_anchor=(1, 1))
    save_current_figure(RESULTS_DIR / f"atac_{patient}_clone_composition.png")


def annotate_patient(patient: str, anno: pd.DataFrame, diploid_cutoff: float = 0.05) -> None:
    path = INTERIM_DIR / f"atac_{patient}.h5ad"
    if not path.exists():
        raise FileNotFoundError(f"Missing ATAC input for {patient}: {path}")

    print(f"[atac] {patient}: reading {path}")
    adata = sc.read_h5ad(path)
    adata.X = adata.X - 1
    adata.obs = adata.obs.join(anno, how="left")
    adata = adata[~adata.obs["cell_type"].isna()].copy()
    adata = add_cnv_metadata(adata)

    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_pcs=min(40, adata.obsm["X_pca"].shape[1]))
    sc.tl.leiden(adata, key_added="cnv_leiden")
    sc.tl.paga(adata, groups="cnv_leiden")
    sc.tl.umap(adata, init_pos="paga")
    cnv.tl.cnv_score(adata)

    clone_scores = adata.obs.groupby("cnv_leiden")["cnv_score"].mean()
    diploid_clones = clone_scores[clone_scores < diploid_cutoff].index
    adata.obs["cnv_leiden"] = adata.obs["cnv_leiden"].replace(
        {clone: "diploid" for clone in diploid_clones}
    )
    adata.uns.pop("dendrogram_cnv_leiden", None)
    adata.uns.pop("cnv_leiden_colors", None)

    plot_qc(adata, patient)

    out = INTERIM_DIR / f"atac_annotated_{patient}.h5ad"
    adata.write_h5ad(out)
    print(f"[atac] {patient}: wrote {out}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--patients", default="all", help="Comma-separated patients or 'all'.")
    parser.add_argument("--diploid-cutoff", type=float, default=0.05)
    args = parser.parse_args()

    anno = load_annotation()
    for patient in parse_patients(args.patients):
        annotate_patient(patient, anno, args.diploid_cutoff)


if __name__ == "__main__":
    main()
