-- mart_signatures.sql
-- =====================
-- Mart: differentially expressed genes between TUMOR and NORMAL samples.
-- Joins to sample metadata to get biological condition labels.
-- Output is the primary table for dashboards and downstream analysis.
--
-- Method: mean log2TPM difference per gene (proxy for log2 fold change).
--         Rank by absolute difference — highest = most changed gene.
--
-- Grain: one row per gene, summarizing expression across all samples.

{{ config(materialized='table', cluster_by=['direction', 'expression_class']) }}

with

expr as (
    select * from {{ ref('int_normalized') }}
),

meta as (
    select
        sample_id,
        -- Standardize condition labels
        case
            when upper(condition) in ('TUMOR','CANCER','PDAC') then 'TUMOR'
            when upper(condition) in ('NORMAL','ADJACENT','HEALTHY') then 'NORMAL'
            else 'UNKNOWN'
        end as condition
    from {{ source('rna', 'sample_metadata') }}
),

-- Mean expression per gene per condition
per_condition as (
    select
        e.gene_id,
        m.condition,
        avg(e.log2_tpm_plus1)   as mean_log2tpm,
        avg(e.tpm)              as mean_tpm,
        count(distinct e.sample_id) as n_samples
    from expr e
    join meta m on e.sample_id = m.sample_id
    where m.condition in ('TUMOR', 'NORMAL')
    group by 1, 2
),

tumor  as (select * from per_condition where condition = 'TUMOR'),
normal as (select * from per_condition where condition = 'NORMAL')

select
    t.gene_id,

    -- Expression in each condition
    t.mean_log2tpm              as mean_log2tpm_tumor,
    n.mean_log2tpm              as mean_log2tpm_normal,
    t.mean_tpm                  as mean_tpm_tumor,
    n.mean_tpm                  as mean_tpm_normal,
    t.n_samples                 as n_tumor_samples,
    n.n_samples                 as n_normal_samples,

    -- Log2 fold change: positive = higher in tumor
    t.mean_log2tpm - n.mean_log2tpm         as log2_fold_change,

    -- Direction
    case
        when (t.mean_log2tpm - n.mean_log2tpm) >  0.5 then 'UP'
        when (t.mean_log2tpm - n.mean_log2tpm) < -0.5 then 'DOWN'
        else 'UNCHANGED'
    end                                     as direction,

    -- Magnitude bucket (useful for volcano plot coloring)
    case
        when abs(t.mean_log2tpm - n.mean_log2tpm) >= 2 then 'STRONG (≥4x)'
        when abs(t.mean_log2tpm - n.mean_log2tpm) >= 1 then 'MODERATE (≥2x)'
        else                                                  'WEAK (<2x)'
    end                                     as fc_bucket,

    -- Expression class in tumor (for filtering)
    case
        when t.mean_tpm >= 100 then 'HIGH'
        when t.mean_tpm >= 10  then 'MEDIUM'
        when t.mean_tpm >= 1   then 'LOW'
        else                        'VERY_LOW'
    end                                     as expression_class,

    -- Rank by magnitude of change (1 = most changed gene)
    row_number() over (
        order by abs(t.mean_log2tpm - n.mean_log2tpm) desc
    )                                       as rank_by_change,

    current_timestamp()                     as dbt_updated_at

from tumor t
join normal n on t.gene_id = n.gene_id
where abs(t.mean_log2tpm - n.mean_log2tpm) > 0.1   -- remove unchanged genes
order by rank_by_change
