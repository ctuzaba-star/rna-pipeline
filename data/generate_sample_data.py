"""
generate_sample_data.py
========================
Generates synthetic RNA-seq data for local development and testing.
No NCBI account or internet connection required.

The synthetic data uses the same schema as real GEO data, with:
  - Realistic expression distributions (log-normal)
  - Known pancreatic cancer genes injected with real fold changes
  - 2 intentionally bad samples for QC testing

Usage:
    python data/generate_sample_data.py

Output:
    data/expression_matrix.parquet  — genes × samples (wide format)
    data/sample_metadata.parquet    — one row per sample
    data/gene_annotation.parquet    — gene symbols and biotypes
"""

import os
from datetime import datetime

import numpy as np
import pandas as pd

rng = np.random.default_rng(seed=42)

# Known pancreatic cancer genes with realistic fold changes
# Source: TCGA PAAD cohort, Bailey et al. 2016 Nature (public data)
KNOWN_GENES = {
    # Upregulated in tumor
    "MUC1":   +4.2,  "MUC5AC": +5.1,  "KRT19":  +4.0,
    "EPCAM":  +3.3,  "MKI67":  +2.5,  "TOP2A":  +2.8,
    "VEGFA":  +2.4,  "MMP9":   +3.1,  "PCNA":   +2.1,
    # Downregulated in tumor (normal pancreas enzymes)
    "PRSS1":  -3.2,  "CELA3A": -4.1,  "PNLIP":  -3.8,
    "AMY2A":  -2.9,  "CPA1":   -3.5,  "INS":    -5.0,
    # Housekeeping (stable — should not change)
    "ACTB":    0.0,  "GAPDH":   0.0,  "B2M":     0.0,
    # Mitochondrial (high = poor sample quality)
    "MT-ND1":  0.2,  "MT-CO1":  0.2,  "MT-CYB":  0.2,
}

N_EXTRA_GENES = 480     # fill out to ~500 total genes
N_TUMOR       = 15
N_NORMAL      = 15


def generate_expression_matrix() -> pd.DataFrame:
    """Generate wide expression matrix: rows=genes, cols=GSM sample IDs."""
    gene_names = list(KNOWN_GENES.keys())
    gene_names += [f"ENSG{str(i).zfill(11)}" for i in range(N_EXTRA_GENES)]
    n_genes = len(gene_names)

    sample_ids = (
        [f"GSM_T{str(i+1).zfill(3)}" for i in range(N_TUMOR)] +
        [f"GSM_N{str(i+1).zfill(3)}" for i in range(N_NORMAL)]
    )
    n_samples = len(sample_ids)

    # Base expression: log-normal (most genes low, few highly expressed)
    base = rng.lognormal(mean=5.0, sigma=2.5, size=(n_genes, n_samples))

    # Apply fold changes for known genes
    for i, gene in enumerate(gene_names):
        if gene in KNOWN_GENES:
            log2fc = KNOWN_GENES[gene]
            # Tumor samples get +fc/2, normal samples get -fc/2
            base[i, :N_TUMOR]  *= 2 ** (log2fc / 2)
            base[i, N_TUMOR:]  *= 2 ** (-log2fc / 2)

    # Degrade 2 samples for QC testing
    mito_indices = [i for i, g in enumerate(gene_names) if g.startswith("MT-")]
    # Sample T001: high mitochondrial fraction (simulate dying cells)
    base[mito_indices, 0] *= 10
    # Sample T002: very low library size (simulate failed sequencing)
    base[:, 1] *= 0.03

    # Round to integers (RNA-seq counts are discrete)
    counts = np.round(base).astype(int)

    df = pd.DataFrame(counts, columns=sample_ids)
    df.insert(0, "gene_id", gene_names)
    return df


def generate_metadata() -> pd.DataFrame:
    """Generate sample metadata with condition labels."""
    rows = []
    for i in range(N_TUMOR):
        rows.append({
            "sample_id":    f"GSM_T{str(i+1).zfill(3)}",
            "condition":    "TUMOR",
            "tissue":       "PANCREAS",
            "organism":     "homo sapiens",
            "platform_id":  "GPL24676",
            "geo_accession": "GSE245552_SYNTHETIC",
        })
    for i in range(N_NORMAL):
        rows.append({
            "sample_id":    f"GSM_N{str(i+1).zfill(3)}",
            "condition":    "NORMAL",
            "tissue":       "PANCREAS",
            "organism":     "homo sapiens",
            "platform_id":  "GPL24676",
            "geo_accession": "GSE245552_SYNTHETIC",
        })
    return pd.DataFrame(rows)


def generate_annotations(gene_ids: list) -> pd.DataFrame:
    """Generate gene annotation table."""
    rows = []
    for gene_id in gene_ids:
        biotype = (
            "mitochondrial" if gene_id.startswith("MT-") else
            "protein_coding" if not gene_id.startswith("ENSG") else
            rng.choice(["protein_coding", "lncRNA", "pseudogene"], p=[0.7, 0.2, 0.1])
        )
        rows.append({
            "gene_id":      gene_id,
            "gene_symbol":  gene_id if not gene_id.startswith("ENSG") else f"GENE_{gene_id[-5:]}",
            "gene_biotype": biotype,
            "chromosome":   rng.choice([str(i) for i in range(1, 23)] + ["X", "Y", "MT"]),
        })
    return pd.DataFrame(rows)


if __name__ == "__main__":
    os.makedirs("data", exist_ok=True)

    print("🧬 Generating synthetic RNA-seq data...")
    print(f"   {N_TUMOR + N_EXTRA_GENES + len(KNOWN_GENES)} genes")
    print(f"   {N_TUMOR} tumor + {N_NORMAL} normal samples")
    print(f"   Based on GSE245552 structure (NCBI GEO, public domain)\n")

    expr_df  = generate_expression_matrix()
    meta_df  = generate_metadata()
    annot_df = generate_annotations(expr_df["gene_id"].tolist())

    expr_df.to_parquet("data/expression_matrix.parquet",  index=False)
    meta_df.to_parquet("data/sample_metadata.parquet",    index=False)
    annot_df.to_parquet("data/gene_annotation.parquet",   index=False)

    print("✅ data/expression_matrix.parquet")
    print("✅ data/sample_metadata.parquet")
    print("✅ data/gene_annotation.parquet")
    print()
    print("Note: GSM_T001 has high mito fraction — should fail QC")
    print("      GSM_T002 has very low library size — should fail QC")
    print()
    print("Next: python spark/normalize.py --input data/expression_matrix.parquet")
