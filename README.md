RNA-seq Pipeline

Processes public RNA-seq data from NCBI GEO through ingestion, normalization, and analysis.

Handles pancreatic cancer vs. normal tissue samples (GSE245552) with TPM normalization and quality control.

---

## Tech Stack

| Component | Technology |
|---|---|
| Ingestion | Python + GEOparse |
| Processing | PySpark (TPM normalization, QC) |
| Warehouse | Snowflake |
| Transforms | dbt |
| Orchestration | Apache Airflow |
| CI/CD | GitHub Actions |

---

## Dataset

**[GSE245552](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245552)** - Pancreatic ductal adenocarcinoma vs. normal tissue
- 62 samples · ~20,000 genes · Illumina NovaSeq · Public domain

---

## Project Structure

```
rna-pipeline/
├── ingestion/
│   └── geo_fetcher.py          # Download data from NCBI GEO
├── spark/
│   └── normalize.py            # TPM normalization + QC in PySpark
├── dbt/
│   └── models/
│       ├── stg_expression.sql  # Clean raw data
│       ├── int_normalized.sql  # TPM + QC flags
│       └── mart_signatures.sql # Top differentially expressed genes
├── data/
│   └── generate_sample_data.py # Generate synthetic test data
├── tests/
│   └── test_normalize.py       # Unit tests
├── .github/workflows/
│   └── ci.yml                  # GitHub Actions CI
├── requirements.txt
└── README.md
```

---

## Quick Start

No external accounts required for local development.

```bash
# Install dependencies
pip install -r requirements.txt

# Generate test data
python data/generate_sample_data.py

# Run normalization
python spark/normalize.py --input data/expression_matrix.parquet

# Run tests
pytest tests/ -v

# Optional: Download real GEO data
python ingestion/geo_fetcher.py --accession GSE245552
```

---

## Pipeline Flow

```
NCBI GEO
    ↓
geo_fetcher.py    → Downloads expression + metadata
    ↓
normalize.py      → TPM normalization, QC
    ↓
Snowflake         → Raw data storage
    ↓
dbt staging       → Data cleaning
    ↓
dbt intermediate  → QC filtering
    ↓
dbt mart         → Gene signatures
```

---

## TPM Normalization

RNA-seq counts vary by sequencing depth. TPM (Transcripts Per Million) normalizes this:

```
TPM_gene = (raw_count / gene_length) / Σ(raw_count / gene_length) × 1,000,000
```

Makes gene expression comparable across samples.

---

## Test Data

Includes synthetic dataset matching real GEO schema for local testing:

```bash
python data/generate_sample_data.py
# Creates:
# - data/expression_matrix.parquet   (genes × samples)
# - data/sample_metadata.parquet     (sample info)
# - data/gene_annotation.parquet     (gene info)
```

Two samples are intentionally poor quality for QC testing:
- `GSM_T001` - High mitochondrial content (cell death)
- `GSM_T002` - Low sequencing depth (failed run)
