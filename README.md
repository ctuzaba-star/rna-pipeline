# 🧬 RNA-seq Expression Pipeline

An end-to-end data engineering pipeline for public RNA-seq data from **NCBI GEO** — the world's largest gene expression repository.

Processes pancreatic cancer vs. normal tissue samples (GSE245552) through ingestion, normalization, transformation, and a Snowflake data warehouse.

---

## Stack

| Layer | Tool |
|---|---|
| Ingestion | Python + GEOparse (NCBI GEO API) |
| Processing | PySpark — TPM normalization, QC |
| Warehouse | Snowflake |
| Transforms | dbt |
| Orchestration | Apache Airflow |
| CI/CD | GitHub Actions |

---

## Dataset

**[GSE245552](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245552)** — Pancreatic ductal adenocarcinoma vs. normal tissue
- 62 samples · ~20,000 genes · Illumina NovaSeq · NIH-funded, public domain

---

## Project Structure

```
rna-pipeline/
├── ingestion/
│   └── geo_fetcher.py          # Download expression data from NCBI GEO
├── spark/
│   └── normalize.py            # TPM normalization + sample QC in PySpark
├── dbt/
│   └── models/
│       ├── stg_expression.sql  # Staging: clean raw data
│       ├── int_normalized.sql  # Intermediate: TPM + QC flags
│       └── mart_signatures.sql # Mart: top differentially expressed genes
├── data/
│   └── generate_sample_data.py # Generate synthetic data (no download needed)
├── tests/
│   └── test_normalize.py       # Unit tests
├── .github/workflows/
│   └── ci.yml                  # GitHub Actions CI
├── requirements.txt
└── README.md
```

---

## Quick Start

**No Snowflake or AWS account needed for local development.**

```bash
# 1. Clone and install
git clone https://github.com/YOUR_USERNAME/rna-pipeline.git
cd rna-pipeline
pip install -r requirements.txt

# 2. Generate synthetic sample data
python data/generate_sample_data.py

# 3. Run normalization
python spark/normalize.py --input data/expression_matrix.parquet

# 4. Run tests
pytest tests/ -v

# 5. (Optional) Download the real GEO dataset
python ingestion/geo_fetcher.py --accession GSE245552
```

---

## Pipeline Overview

```
NCBI GEO (public)
      │
      ▼
geo_fetcher.py          ← downloads expression matrix + sample metadata
      │
      ▼
normalize.py            ← TPM normalization, log2 transform, sample QC
      │
      ▼
Snowflake RAW_DB        ← COPY INTO from S3/local parquet
      │
      ▼
dbt staging             ← type cast, deduplicate, null guards
      │
      ▼
dbt intermediate        ← QC flags, filter low-expressed genes
      │
      ▼
dbt mart                ← top differentially expressed genes, ranked
```

---

## What is TPM Normalization?

Raw RNA-seq counts are confounded by how many total reads were sequenced per sample. TPM (Transcripts Per Million) fixes this:

```
TPM_gene = (raw_count / gene_length) / Σ(raw_count / gene_length) × 1,000,000
```

This makes gene expression comparable across samples regardless of sequencing depth.

---

## Sample Data

The repo ships with a synthetic dataset (same schema as real GEO data) so you can run everything locally without downloading anything:

```bash
python data/generate_sample_data.py
# → data/expression_matrix.parquet   (genes × samples)
# → data/sample_metadata.parquet     (sample conditions)
# → data/gene_annotation.parquet     (gene symbols, biotypes)
```

Two samples are intentionally degraded to test the QC step:
- `GSM_T001` — high mitochondrial fraction (dying cells)
- `GSM_T002` — very low library size (failed sequencing run)
