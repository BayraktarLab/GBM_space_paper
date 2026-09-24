"""Convert epiAneufinder sample outputs into patient-level ATAC h5ad inputs.

This creates the unannotated files consumed by 01_atac_clone_annotation.py:

    ../../data/interim/atac_<PATIENT>.h5ad

The converter expects epiAneufinder's per-sample `results_table.tsv`, whose
values are copy-number states per cell/bin: 0 loss, 1 normal, 2 gain. The
original ATAC annotation script subtracts 1 from `adata.X`, so this converter
stores the raw 0/1/2 states and lets the downstream script transform them to
-1/0/1.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from pipeline_config import INTERIM_DIR, PROJECT_ROOT, parse_patients


DEFAULT_EPI_ROOT = PROJECT_ROOT / "atac_cn_calling"


def find_results_table(sample_outdir: Path) -> Path:
    candidates = [
        sample_outdir / "results_table.tsv",
        sample_outdir / "epiAneufinder_results" / "results_table.tsv",
    ]
    candidates.extend(sample_outdir.rglob("results_table.tsv"))
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f"No results_table.tsv found under {sample_outdir}")


def parse_region(region: str) -> tuple[str, int, int]:
    """Parse common epiAneufinder/bin labels into chromosome/start/end."""
    text = str(region)
    patterns = [
        r"^(?P<chrom>[^:_-]+):(?P<start>[0-9]+)-(?P<end>[0-9]+)$",
        r"^(?P<chrom>chr[^:_-]+)_(?P<start>[0-9]+)_(?P<end>[0-9]+)$",
        r"^(?P<chrom>chr[^:_-]+)-(?P<start>[0-9]+)-(?P<end>[0-9]+)$",
    ]
    for pattern in patterns:
        match = re.match(pattern, text)
        if match:
            return (
                match.group("chrom"),
                int(match.group("start")),
                int(match.group("end")),
            )
    raise ValueError(f"Could not parse genomic bin label: {region}")


def normalize_results_table(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep=r"\s+", index_col=0)
    table = table.dropna(how="all")

    lower_columns = {str(col).lower(): col for col in table.columns}
    chrom_col = (
        lower_columns.get("seq")
        or lower_columns.get("chr")
        or lower_columns.get("chromosome")
        or lower_columns.get("seqnames")
    )
    start_col = lower_columns.get("start")
    end_col = lower_columns.get("end")
    if chrom_col is not None and start_col is not None and end_col is not None:
        regions = (
            table[chrom_col].astype(str)
            + ":"
            + table[start_col].astype(int).astype(str)
            + "-"
            + table[end_col].astype(int).astype(str)
        )
        value_table = table.drop(columns=[chrom_col, start_col, end_col])
        value_table.index = regions
        table = value_table.T

    first_index = str(table.index[0]) if len(table.index) else ""
    index_looks_like_bins = first_index.startswith("chr") and any(token in first_index for token in [":", "_", "-"])
    column_looks_like_bins = any(str(col).startswith("chr") for col in table.columns[: min(5, len(table.columns))])

    if index_looks_like_bins and not column_looks_like_bins:
        table = table.T

    # Keep numeric copy-number states only. This avoids metadata columns if the
    # result table contains extras.
    numeric = table.apply(pd.to_numeric, errors="coerce")
    numeric = numeric.dropna(axis=1, how="all")
    numeric = numeric.dropna(axis=0, how="all")
    if numeric.empty:
        raise ValueError(f"No numeric copy-number matrix found in {path}")
    return numeric


def table_to_anndata(table: pd.DataFrame, patient: str, section: str, sample: str) -> ad.AnnData:
    bins = list(table.columns)
    parsed_bins = [parse_region(region) for region in bins]
    var = pd.DataFrame(parsed_bins, columns=["chromosome", "start", "end"], index=bins)
    var = var.sort_values(["chromosome", "start", "end"], kind="mergesort")
    table = table.loc[:, var.index]

    barcodes = table.index.astype(str).str.replace(r"^cell-", "", regex=True)
    obs = pd.DataFrame(index=[f"{sample}#{barcode}" for barcode in barcodes])
    obs["barcode"] = barcodes
    obs["patient"] = patient
    obs["section"] = section
    obs["sample"] = sample

    matrix = sparse.csr_matrix(table.to_numpy(dtype=np.float32))
    return ad.AnnData(X=matrix, obs=obs, var=var)


def find_epianeufinder_outputs(epi_root: Path, resolution: str) -> pd.DataFrame:
    tables = sorted(epi_root.glob(f"AT*-*/cellranger-arc*/{resolution}/epiAneufinder_results/results_table.tsv"))
    if not tables:
        tables = sorted(epi_root.glob(f"AT*-*/*/{resolution}/epiAneufinder_results/results_table.tsv"))
    if not tables:
        raise FileNotFoundError(
            f"No {resolution}/epiAneufinder_results/results_table.tsv files found under {epi_root}"
        )

    rows = []
    for result_path in tables:
        section = result_path.parents[3].name
        sample = result_path.parents[2].name
        patient = section.split("-")[0]
        rows.append(
            {
                "patient": patient,
                "section": section,
                "sample": sample,
                "epianeufinder_outdir": str(result_path.parent),
            }
        )
    return pd.DataFrame(rows)


def build_patient_atac(patient: str, outputs: pd.DataFrame, force: bool) -> None:
    patient_rows = outputs[outputs["patient"].astype(str).str.upper() == patient]
    if patient_rows.empty:
        print(f"[convert] {patient}: no epiAneufinder outputs found")
        return

    out_path = INTERIM_DIR / f"atac_{patient}.h5ad"
    if out_path.exists() and not force:
        print(f"[convert] {patient}: {out_path} exists; use --force to overwrite")
        return

    adatas = []
    for row in patient_rows.itertuples(index=False):
        result_path = find_results_table(Path(row.epianeufinder_outdir))
        print(f"[convert] {patient}: reading {result_path}")
        table = normalize_results_table(result_path)
        adatas.append(table_to_anndata(table, patient, row.section, row.sample))

    combined = ad.concat(adatas, join="outer", merge="same")
    combined.var = combined.var.astype({"start": int, "end": int})
    combined.write_h5ad(out_path)
    print(f"[convert] {patient}: wrote {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--epi-root", type=Path, default=DEFAULT_EPI_ROOT)
    parser.add_argument("--resolution", default="1mb", help="epiAneufinder bin-size folder to read.")
    parser.add_argument("--patients", default="all", help="Comma-separated patients or 'all'.")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    epi_root = args.epi_root
    if not epi_root.is_absolute():
        epi_root = (PROJECT_ROOT / epi_root).resolve() if str(epi_root).startswith("data/") else epi_root.resolve()

    outputs = find_epianeufinder_outputs(epi_root, args.resolution)
    for patient in parse_patients(args.patients):
        build_patient_atac(patient, outputs, args.force)


if __name__ == "__main__":
    main()
