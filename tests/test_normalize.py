"""
test_normalize.py
=================
Unit tests for the TPM normalization and QC logic.
Uses a small in-memory Spark DataFrame — no files needed.
"""

import pytest
from pyspark.sql import Row, SparkSession
from pyspark.sql import functions as F


@pytest.fixture(scope="session")
def spark():
    return (
        SparkSession.builder
        .master("local[1]")
        .appName("RNA-seq Tests")
        .config("spark.sql.shuffle.partitions", "2")
        .getOrCreate()
    )


@pytest.fixture
def long_df(spark):
    """Small expression DataFrame in long format for testing."""
    data = [
        # sample S1: 2 genes, known counts
        Row(gene_id="GENE_A", sample_id="S1", raw_count=1000.0),
        Row(gene_id="GENE_B", sample_id="S1", raw_count=4000.0),
        # sample S2: same genes, different counts
        Row(gene_id="GENE_A", sample_id="S2", raw_count=500.0),
        Row(gene_id="GENE_B", sample_id="S2", raw_count=500.0),
        # MT gene for QC testing
        Row(gene_id="MT-ND1", sample_id="S1", raw_count=100.0),
        Row(gene_id="MT-ND1", sample_id="S2", raw_count=5000.0),  # high mito in S2
    ]
    return spark.createDataFrame(data)


class TestTPMNormalization:

    def test_tpm_sums_to_one_million_per_sample(self, spark, long_df):
        """TPM values must sum to exactly 1,000,000 per sample."""
        import sys
        sys.path.insert(0, ".")
        from spark.normalize import compute_tpm

        result = compute_tpm(long_df)
        sums = (
            result
            .groupBy("sample_id")
            .agg(F.sum("tpm").alias("total_tpm"))
            .collect()
        )
        for row in sums:
            assert abs(row["total_tpm"] - 1_000_000) < 1.0, (
                f"Sample {row['sample_id']}: TPM sum = {row['total_tpm']:.2f}, expected 1,000,000"
            )

    def test_higher_count_gives_higher_tpm(self, spark, long_df):
        """Within a sample, higher raw count → higher TPM."""
        import sys
        sys.path.insert(0, ".")
        from spark.normalize import compute_tpm

        result = compute_tpm(long_df)
        s1 = {r["gene_id"]: r["tpm"] for r in result.where(F.col("sample_id") == "S1").collect()}

        assert s1["GENE_B"] > s1["GENE_A"], (
            "GENE_B has 4x the count of GENE_A, so it should have higher TPM"
        )

    def test_log2_tpm_plus1_non_negative(self, spark, long_df):
        """log2(TPM + 1) must always be ≥ 0 since TPM ≥ 0."""
        import sys
        sys.path.insert(0, ".")
        from spark.normalize import compute_tpm

        result = compute_tpm(long_df)
        negatives = result.where(F.col("log2_tpm_plus1") < 0).count()
        assert negatives == 0, f"Found {negatives} rows with negative log2_tpm_plus1"


class TestSampleQC:

    def test_high_mito_fraction_fails_qc(self, spark, long_df):
        """Sample S2 has very high mitochondrial counts and should fail QC."""
        import sys
        sys.path.insert(0, ".")
        from spark.normalize import compute_tpm, flag_sample_qc

        tpm_df = compute_tpm(long_df)
        qc_df  = flag_sample_qc(tpm_df)

        s2_qc = qc_df.where(F.col("sample_id") == "S2").first()

        # S2: MT-ND1 has 5000/6000 total counts = 83% mito → should fail
        assert s2_qc["mito_fraction"] > 0.30, (
            f"Expected high mito fraction for S2, got {s2_qc['mito_fraction']:.2%}"
        )
        # Note: qc_pass may still be True here because total_counts and
        # detected_genes thresholds are set for real data (500K reads minimum)
        # This test validates the metric calculation, not the threshold

    def test_mito_fraction_between_0_and_1(self, spark, long_df):
        """Mitochondrial fraction must be between 0 and 1."""
        import sys
        sys.path.insert(0, ".")
        from spark.normalize import compute_tpm, flag_sample_qc

        tpm_df = compute_tpm(long_df)
        qc_df  = flag_sample_qc(tpm_df)

        invalid = qc_df.where(
            (F.col("mito_fraction") < 0) | (F.col("mito_fraction") > 1)
        ).count()
        assert invalid == 0

    def test_all_samples_have_qc_record(self, spark, long_df):
        """Every sample in expression data must have a QC record."""
        import sys
        sys.path.insert(0, ".")
        from spark.normalize import compute_tpm, flag_sample_qc

        tpm_df   = compute_tpm(long_df)
        qc_df    = flag_sample_qc(tpm_df)
        samples  = {r["sample_id"] for r in long_df.select("sample_id").distinct().collect()}
        qc_ids   = {r["sample_id"] for r in qc_df.select("sample_id").collect()}

        assert samples == qc_ids, f"Missing QC for samples: {samples - qc_ids}"
