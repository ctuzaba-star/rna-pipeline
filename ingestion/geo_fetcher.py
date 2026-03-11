"""
geo_fetcher.py
==============
Downloads RNA-seq expression data from NCBI GEO (Gene Expression Omnibus).

Usage:
    python ingestion/geo_fetcher.py --accession GSE245552

What it does:
    1. Downloads the series matrix file from NCBI FTP
    2. Parses gene expression values and sample metadata
    3. Saves as Parquet files ready for the normalize.py step

Dataset: GSE245552 — Pancreatic cancer vs normal tissue (public domain)
"""

import argparse
import gzip
import logging
import os
from pathlib import Path

import pandas as pd
import requests

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s")
log = logging.getLogger(__name__)

NCBI_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def download_matrix(accession: str, output_dir: str) -> Path:
    """Download the GSE series matrix file from NCBI FTP."""
    # GEO FTP path uses a 'nnn' prefix: GSE245552 → GSE245nnn
    prefix = accession[:-3] + "nnn"
    url = f"{NCBI_FTP}/{prefix}/{accession}/matrix/{accession}_series_matrix.txt.gz"

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"{accession}_series_matrix.txt.gz"

    if dest.exists():
        log.info("Already downloaded: %s", dest.name)
        return dest

    log.info("Downloading %s from NCBI...", accession)
    log.info("URL: %s", url)

    resp = requests.get(url, stream=True, timeout=120)
    resp.raise_for_status()

    total_mb = 0
    with open(dest, "wb") as f:
        for chunk in resp.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)
            total_mb += len(chunk) / 1e6

    log.info("Downloaded %.1f MB → %s", total_mb, dest)
    return dest


def parse_matrix(matrix_gz: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse SOFT matrix file into expression and metadata DataFrames.

    The file has two sections:
      - Header lines starting with '!' = sample metadata
      - Table between !series_matrix_table_begin / end = expression values
    """
    sample_meta = {}
    in_table = False
    header = None
    data_rows = []

    opener = gzip.open if str(matrix_gz).endswith(".gz") else open

    with opener(matrix_gz, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")

            # Collect sample metadata
            if line.startswith("!Sample_geo_accession"):
                for gsm in line.split("\t")[1:]:
                    gsm = gsm.strip('"')
                    sample_meta[gsm] = {"sample_id": gsm}

            elif line.startswith("!Sample_") and sample_meta:
                key = line.split("=")[0].lstrip("!").strip()
                vals = line.split("\t")[1:]
                for i, gsm in enumerate(list(sample_meta.keys())):
                    if i < len(vals):
                        sample_meta[gsm][key] = vals[i].strip('" ')

            # Parse expression table
            elif line == "!series_matrix_table_begin":
                in_table = True

            elif line == "!series_matrix_table_end":
                in_table = False

            elif in_table:
                if header is None:
                    header = [h.strip('"') for h in line.split("\t")]
                else:
                    data_rows.append(line.split("\t"))

    # Build expression DataFrame
    if header and data_rows:
        expr_df = pd.DataFrame(
            [r for r in data_rows if len(r) == len(header)],
            columns=header,
        )
        expr_df = expr_df.rename(columns={header[0]: "gene_id"})
        gsm_cols = [c for c in expr_df.columns if c.startswith("GSM")]
        for col in gsm_cols:
            expr_df[col] = pd.to_numeric(expr_df[col], errors="coerce")
    else:
        expr_df = pd.DataFrame(columns=["gene_id"])

    # Build metadata DataFrame
    meta_df = pd.DataFrame(list(sample_meta.values()))
    if not meta_df.empty:
        meta_df["condition"] = meta_df.apply(_infer_condition, axis=1)

    return expr_df, meta_df


def _infer_condition(row) -> str:
    """Infer tumor/normal from free-text GEO metadata fields."""
    text = " ".join(str(v) for v in row.values()).lower()
    if any(w in text for w in ["tumor", "cancer", "pdac", "malignant"]):
        return "TUMOR"
    elif any(w in text for w in ["normal", "adjacent", "healthy"]):
        return "NORMAL"
    return "UNKNOWN"


def fetch(accession: str, output_dir: str = "data/raw") -> None:
    """Main entry point: download and save GEO dataset as Parquet."""
    matrix_gz = download_matrix(accession, output_dir)

    log.info("Parsing expression matrix...")
    expr_df, meta_df = parse_matrix(matrix_gz)

    n_genes   = len(expr_df)
    n_samples = len([c for c in expr_df.columns if c.startswith("GSM")])
    log.info("Parsed: %d genes × %d samples", n_genes, n_samples)

    # Save as Parquet
    expr_path = Path(output_dir) / "expression_matrix.parquet"
    meta_path = Path(output_dir) / "sample_metadata.parquet"

    expr_df.to_parquet(expr_path, index=False)
    meta_df.to_parquet(meta_path, index=False)

    log.info("Saved expression matrix → %s", expr_path)
    log.info("Saved sample metadata  → %s", meta_path)
    log.info("Done! Next step: python spark/normalize.py --input %s", output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--accession", default="GSE245552")
    parser.add_argument("--output",    default="data/raw")
    args = parser.parse_args()
    fetch(args.accession, args.output)
