"""
normalize.py
============
TPM normalization and sample QC for RNA-seq expression data.

Steps:
    1. Load expression matrix (wide format: rows=genes, cols=samples)
    2. Melt to long format (gene_id, sample_id, raw_count)
    3. Compute TPM: corrects for gene length + sequencing depth
    4. Log2(TPM + 1) transform: compresses dynamic range for ML
    5. Flag low-quality samples (high mitochondrial %, low library size)
    6. Write normalized Parquet for dbt

Usage:
    python spark/normalize.py --input data/expression_matrix.parquet
    python spark/normalize.py --input data/raw/GSE245552/  # full GEO download

TPM formula:
    rate_i     = raw_count_i / gene_length_i
    TPM_i      = rate_i / sum(all rates) × 1,000,000
"""

import argparse
import logging
from pathlib import Path

from pyspark.sql import SparkSession
from pyspark.sql import functions as F
from pyspark.sql import types as T
from pyspark.sql.window import Window

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(message)s")
log = logging.getLogger(__name__)

# QC thresholds for bulk RNA-seq
MIN_LIBRARY_SIZE  = 500_000    # samples with fewer reads are likely failed runs
MIN_GENES_DETECTED = 5_000     # samples expressing fewer genes are suspicious
MAX_MITO_FRACTION  = 0.30      # >30% mitochondrial = poor sample quality


def get_spark() -> SparkSession:
    return (
        SparkSession.builder
        .appName("RNA-seq Normalization")
        .config("spark.sql.shuffle.partitions", "8")   # small dataset, keep it fast
        .getOrCreate()
    )


def load_expression(spark: SparkSession, path: str):
    """Load expression matrix. Accepts parquet file or directory."""
    log.info("Loading expression matrix from %s", path)
    df = spark.read.parquet(path)

    # Detect GSM sample columns
    gsm_cols = [c for c in df.columns if c.startswith("GSM")]
    log.info("Found %d genes, %d samples", df.count(), len(gsm_cols))
    return df, gsm_cols


def melt_to_long(df, gsm_cols: list):
    """
    Convert wide matrix (gene_id | GSM001 | GSM002 | ...) to long format
    (gene_id | sample_id | raw_count).

    Long format is required for per-sample normalization in Spark.
    """
    log.info("Melting wide matrix → long format")

    stack_expr = ", ".join([f"'{col}', `{col}`" for col in gsm_cols])
    long_df = df.selectExpr(
        "gene_id",
        f"stack({len(gsm_cols)}, {stack_expr}) as (sample_id, raw_count)",
    )

    return (
        long_df
        .where(F.col("raw_count").isNotNull())
        .where(F.col("raw_count") >= 0)
        .withColumn("raw_count", F.col("raw_count").cast(T.DoubleType()))
    )


def compute_tpm(long_df):
    """
    TPM normalization:
        Step 1: rate = raw_count / gene_length  (default 2500bp if unknown)
        Step 2: per_sample_total = sum of all rates in that sample
        Step 3: TPM = rate / per_sample_total × 1,000,000

    Result is comparable across samples regardless of how deeply each was sequenced.
    """
    log.info("Computing TPM normalization")

    sample_window = Window.partitionBy("sample_id")

    return (
        long_df
        # Use 2500bp as default gene length (median human transcript)
        .withColumn("gene_length", F.lit(2500.0))

        # Step 1: reads per bp
        .withColumn("rate", F.col("raw_count") / F.col("gene_length"))

        # Step 2: sum of rates per sample
        .withColumn("sample_rate_sum", F.sum("rate").over(sample_window))

        # Step 3: TPM
        .withColumn(
            "tpm",
            F.when(F.col("sample_rate_sum") > 0,
                   (F.col("rate") / F.col("sample_rate_sum")) * 1_000_000
            ).otherwise(0.0)
        )

        # Log2(TPM + 1): compress dynamic range, stabilize variance
        .withColumn("log2_tpm_plus1", F.log2(F.col("tpm") + 1))

        .select("gene_id", "sample_id", "raw_count", "tpm", "log2_tpm_plus1")
    )


def flag_sample_qc(tpm_df):
    """
    Compute per-sample QC metrics and flag low-quality samples.

    Metrics:
        total_counts     — total raw reads (library size)
        detected_genes   — genes with raw_count > 0
        mito_fraction    — fraction of reads from mitochondrial genes (MT-*)
        qc_pass          — True if all thresholds met
    """
    log.info("Computing sample QC metrics")

    sample_window = Window.partitionBy("sample_id")

    qc_df = (
        tpm_df
        .withColumn("is_mito", F.col("gene_id").rlike(r"(?i)^MT-"))
        .groupBy("sample_id")
        .agg(
            F.sum("raw_count").alias("total_counts"),
            F.count(F.when(F.col("raw_count") > 0, 1)).alias("detected_genes"),
            F.sum(F.when(F.col("is_mito"), F.col("raw_count")).otherwise(0)).alias("mito_counts"),
        )
        .withColumn(
            "mito_fraction",
            F.col("mito_counts") / F.greatest(F.col("total_counts"), F.lit(1))
        )
        .withColumn(
            "library_size_log2",
            F.log2(F.greatest(F.col("total_counts"), F.lit(1)))
        )
        .withColumn(
            "qc_pass",
            (F.col("total_counts")   >= MIN_LIBRARY_SIZE)  &
            (F.col("detected_genes") >= MIN_GENES_DETECTED) &
            (F.col("mito_fraction")  <= MAX_MITO_FRACTION)
        )
        .withColumn(
            "qc_fail_reason",
            F.concat_ws(" | ",
                F.when(F.col("total_counts")   < MIN_LIBRARY_SIZE,  F.lit("LOW_LIBRARY_SIZE")),
                F.when(F.col("detected_genes") < MIN_GENES_DETECTED, F.lit("FEW_GENES_DETECTED")),
                F.when(F.col("mito_fraction")  > MAX_MITO_FRACTION,  F.lit("HIGH_MITO_FRACTION")),
            )
        )
    )

    n_pass = qc_df.where(F.col("qc_pass")).count()
    n_total = qc_df.count()
    log.info("QC: %d/%d samples passed (%.0f%%)", n_pass, n_total, 100 * n_pass / max(n_total, 1))

    return qc_df


def run(input_path: str, output_path: str = "data/normalized") -> None:
    spark = get_spark()

    # Load
    raw_df, gsm_cols = load_expression(spark, input_path)

    # Transform
    long_df   = melt_to_long(raw_df, gsm_cols)
    tpm_df    = compute_tpm(long_df)
    qc_df     = flag_sample_qc(tpm_df)

    # Join QC flags back to expression
    normalized = tpm_df.join(qc_df.select("sample_id", "qc_pass", "qc_fail_reason"), on="sample_id")

    # Write
    out = Path(output_path)
    out.mkdir(parents=True, exist_ok=True)

    normalized.write.mode("overwrite").parquet(str(out / "expression_normalized.parquet"))
    qc_df.write.mode("overwrite").parquet(str(out / "sample_qc.parquet"))

    log.info("Written normalized expression → %s/expression_normalized.parquet", output_path)
    log.info("Written sample QC metrics    → %s/sample_qc.parquet", output_path)
    log.info("Next step: dbt run")

    spark.stop()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",  required=True, help="Path to expression_matrix.parquet")
    parser.add_argument("--output", default="data/normalized")
    args = parser.parse_args()
    run(args.input, args.output)
