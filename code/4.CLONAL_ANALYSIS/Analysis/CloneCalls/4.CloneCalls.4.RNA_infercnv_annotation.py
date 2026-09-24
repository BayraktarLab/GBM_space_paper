"""Run RNA inferCNV clone annotation for all GBM patients.

Inputs
------
data/raw/multiome/*/*/filtered_feature_bc_matrix.h5
data/interim/annotation_obs.29_06_23.csv
data/interim/gene_list

Outputs
-------
data/interim/<patient_lower>_annotation_5k.h5ad
notebooks/results/rna_<PATIENT>_*.png
"""

from __future__ import annotations

import argparse
import glob
import warnings

import anndata as ad
import infercnvpy as cnv
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm

from pipeline_config import INTERIM_DIR, RAW_DIR, REFERENCE_CELL_TYPES, RESULTS_DIR, parse_patients


warnings.simplefilter("ignore")
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def save_current_figure(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight", dpi=180)
    plt.close("all")


def load_gene_positions() -> pd.DataFrame:
    gene_list = pd.read_csv(INTERIM_DIR / "gene_list", sep="\t", header=None)
    genes = []
    for tup in gene_list.itertuples(index=False):
        gene = ""
        try:
            info = tup[8].split(";")
            for item in info:
                label, value = item.split("=")
                if label == "gene_name":
                    gene = value
        except (AttributeError, IndexError, ValueError):
            pass
        genes.append(gene)

    gene_list.index = genes
    gene_list = gene_list[[0, 3, 4]]
    gene_list.columns = ["chromosome", "start", "end"]
    gene_list = gene_list.drop_duplicates()
    return gene_list.reset_index().drop_duplicates(subset="index", keep="last").set_index("index")


def load_annotation() -> pd.DataFrame:
    anno = pd.read_csv(INTERIM_DIR / "annotation_obs.29_06_23.csv")
    anno["cell_type"] = anno["TME_coarse_GBM_granular"].astype(str).str.split(".").str[0]
    return anno


def parse_multiome_path(path: str) -> tuple[str, str]:
    parts = path.split("/")
    section, sample = parts[-3], parts[-2]
    return section, sample


def load_multiome_rna() -> sc.AnnData:
    rna_paths = sorted(glob.glob(str(RAW_DIR / "multiome" / "*" / "*" / "filtered_feature_bc_matrix.h5")))
    if not rna_paths:
        raise FileNotFoundError("No multiome RNA matrices found under data/raw/multiome/*/*/")

    prefix = "/lustre/scratch126/cellgen/team283/gd11/gd11/data_GBM/cellbender_input/data/"
    suffix = "/cellbender_out/_"
    adatas = []
    for path in tqdm(rna_paths, desc="Reading multiome RNA"):
        section, sample = parse_multiome_path(path)
        adata = sc.read_10x_h5(path)
        adata.obs["section"] = section
        adata.obs["sample"] = sample
        adata.obs["cell_id"] = [f"{sample}#{barcode}" for barcode in adata.obs.index]
        adata.obs.index = [f"{prefix}{sample}{suffix}{barcode}" for barcode in adata.obs.index]
        adata.var_names_make_unique()
        adatas.append(adata)

    adata = sc.concat(adatas, axis=0, join="outer")
    adata.var_names_make_unique()
    return adata


def prepare_rna() -> sc.AnnData:
    adata = load_multiome_rna()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata = adata[
        (adata.obs.n_genes_by_counts > 500)
        & (adata.obs.total_counts > 5000)
        & (adata.obs.pct_counts_mt < 10),
        :,
    ].copy()
    adata.raw = adata
    adata.layers["counts"] = adata.X.copy()
    adata.obs["patient"] = adata.obs["section"].astype(str).str.split("-").str[0]

    anno = load_annotation()
    adata.obs = (
        adata.obs.reset_index()
        .merge(anno, on=["cell_id"], how="left")
        .set_index("index")
    )
    adata = adata[~adata.obs["cell_type"].isna()].copy()
    adata = adata[adata.obs["patient"] != "download"].copy()

    gene_positions = load_gene_positions()
    common_genes = adata.var_names.intersection(gene_positions.index)
    adata = adata[:, common_genes].copy()
    adata.var = adata.var.join(gene_positions, how="left")
    return adata


def plot_patient_qc(adata: sc.AnnData, adata_m: sc.AnnData, patient: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    cnv.pl.umap(adata, color="cnv_score", ax=axes[0, 0], show=False)
    cnv.pl.umap(adata, color="cell_type", ax=axes[0, 1], show=False)
    cnv.pl.umap(adata, color="is_reference", ax=axes[1, 0], show=False)
    cnv.pl.umap(adata, color="section", ax=axes[1, 1], show=False)
    save_current_figure(RESULTS_DIR / f"rna_{patient}_infercnv_umap.png")

    if adata_m.n_obs:
        cnv.pl.chromosome_heatmap_summary(
            adata_m, groupby="cnv_leiden", dendrogram=True, figsize=(16, 10), show=False
        )
        save_current_figure(RESULTS_DIR / f"rna_{patient}_malignant_cnv_summary.png")

        df = adata_m.obs[["cell_type", "cnv_leiden", "cnv_score"]].copy()
        df.columns = ["cell_type", "clone", "cnv_score"]
        grouped = df.groupby("clone")["cell_type"].value_counts(normalize=True).unstack().fillna(0)
        cnv_score = df.groupby("clone")["cnv_score"].mean().sort_values(ascending=False)
        ax = grouped.loc[cnv_score.index].plot(kind="bar", stacked=True, figsize=(10, 7))
        ax.plot(range(len(cnv_score)), cnv_score.values * 10, label="10 x cnv_score", ls="--", color="black")
        ax.set_ylabel("Fraction of cells / scaled CNV score")
        ax.legend(bbox_to_anchor=(1, 1))
        save_current_figure(RESULTS_DIR / f"rna_{patient}_clone_composition.png")


def annotate_patient(adata: sc.AnnData, patient: str, malignant_cutoff: float, n_jobs: int) -> None:
    print(f"[rna] {patient}: starting")
    adata_tmp = adata[adata.obs.patient == patient].copy()
    if adata_tmp.n_obs == 0:
        raise ValueError(f"No RNA cells found for {patient}")

    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)
    sc.pp.highly_variable_genes(adata_tmp, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.tl.pca(adata_tmp, svd_solver="arpack")
    sc.pp.neighbors(adata_tmp)
    sc.tl.leiden(adata_tmp)
    sc.tl.paga(adata_tmp)
    sc.tl.umap(adata_tmp, init_pos="paga")

    autosomal = ~adata_tmp.var["chromosome"].isin(["chrM", "chrX", "chrY"])
    adata_tmp = adata_tmp[:, adata_tmp.var[autosomal].index].copy()
    ref_cells = sorted(set(REFERENCE_CELL_TYPES).intersection(adata_tmp.obs.cell_type.unique()))
    if not ref_cells:
        raise ValueError(f"No reference cell types available for {patient}")

    cnv.tl.infercnv(
        adata=adata_tmp,
        reference_key="cell_type",
        reference_cat=ref_cells,
        window_size=250,
        step=1,
        n_jobs=n_jobs,
    )
    cnv.tl.pca(adata_tmp)
    cnv.pp.neighbors(adata_tmp)
    cnv.tl.leiden(adata_tmp)
    cnv.tl.umap(adata_tmp)
    cnv.tl.cnv_score(adata_tmp)
    adata_tmp.obs["is_reference"] = adata_tmp.obs.cell_type.isin(REFERENCE_CELL_TYPES).astype("category")
    adata_tmp.obs["malignant"] = (
        ~adata_tmp.obs.cell_type.isin(REFERENCE_CELL_TYPES)
    ) & (adata_tmp.obs.cnv_score > malignant_cutoff)

    adata_m = adata_tmp[adata_tmp.obs["malignant"]].copy()
    if adata_m.n_obs:
        cnv.tl.pca(adata_m)
        cnv.pp.neighbors(adata_m)
        cnv.tl.leiden(adata_m)
        cnv.tl.umap(adata_m)
        adata_m.uns.pop("dendrogram_cnv_leiden", None)
        adata_m.uns.pop("cnv_leiden_colors", None)

    adata_ref = adata_tmp[~adata_tmp.obs["malignant"]].copy()
    adata_ref.obs["cnv_leiden"] = "diploid"
    to_save = ad.concat([adata_m, adata_ref])
    to_save.uns["cnv"] = adata_tmp.uns["cnv"]
    plot_patient_qc(adata_tmp, adata_m, patient)

    out = INTERIM_DIR / f"{patient.lower()}_annotation_5k.h5ad"
    to_save.write_h5ad(out)
    print(f"[rna] {patient}: wrote {out}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--patients", default="all", help="Comma-separated patients or 'all'.")
    parser.add_argument("--malignant-cutoff", type=float, default=0.004)
    parser.add_argument("--n-jobs", type=int, default=10)
    args = parser.parse_args()

    adata = prepare_rna()
    for patient in parse_patients(args.patients):
        annotate_patient(adata, patient, args.malignant_cutoff, args.n_jobs)


if __name__ == "__main__":
    main()
