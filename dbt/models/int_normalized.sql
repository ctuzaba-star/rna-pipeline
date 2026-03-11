-- int_normalized.sql
-- ====================
-- Intermediate: QC-filtered expression with per-gene summary stats.
-- Only includes samples that passed QC.
-- Adds expression class (HIGH/MEDIUM/LOW) for dashboard filtering.

{{ config(materialized='table') }}

with

passing_samples as (
    -- Only keep samples that passed quality control
    select * from {{ ref('stg_expression') }}
    where qc_pass = true
),

gene_stats as (
    -- Per-gene statistics across all passing samples
    select
        gene_id,
        count(distinct sample_id)               as n_samples,
        avg(tpm)                                as mean_tpm,
        stddev(tpm)                             as std_tpm,
        avg(log2_tpm_plus1)                     as mean_log2tpm,
        sum(case when tpm = 0 then 1 else 0 end) as n_zero_samples
    from passing_samples
    group by gene_id
)

select
    e.gene_id,
    e.sample_id,
    e.raw_count,
    e.tpm,
    e.log2_tpm_plus1,

    -- Expression class (useful for filtering in dashboards)
    case
        when gs.mean_tpm >= 100 then 'HIGH'
        when gs.mean_tpm >= 10  then 'MEDIUM'
        when gs.mean_tpm >= 1   then 'LOW'
        else                         'VERY_LOW'
    end                             as expression_class,

    -- Coefficient of variation: how variable is this gene across samples?
    case
        when gs.mean_tpm > 0
        then gs.std_tpm / gs.mean_tpm
    end                             as coeff_of_variation,

    -- Fraction of samples where gene is not expressed at all
    gs.n_zero_samples::float
        / nullif(gs.n_samples, 0)   as pct_zero_samples,

    -- Mitochondrial gene flag (high expression = cell stress)
    gene_id rlike '(?i)^MT-'        as is_mitochondrial,

    e.dbt_updated_at

from passing_samples e
join gene_stats gs on e.gene_id = gs.gene_id

-- Filter out genes with zero expression in >80% of samples (likely noise)
where gs.n_zero_samples::float / nullif(gs.n_samples, 0) < 0.80
