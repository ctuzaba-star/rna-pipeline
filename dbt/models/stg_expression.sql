-- stg_expression.sql
-- ====================
-- Staging: clean and type-cast the raw normalized expression data.
-- One row per (gene, sample). Removes nulls and negative values.

{{ config(materialized='incremental', unique_key=['gene_id', 'sample_id']) }}

with source as (
    select * from {{ source('rna', 'expression_normalized') }}
    {% if is_incremental() %}
        where loaded_at > (select max(loaded_at) from {{ this }})
    {% endif %}
)

select
    trim(upper(gene_id))                    as gene_id,
    trim(upper(sample_id))                  as sample_id,
    cast(raw_count as double)               as raw_count,
    cast(tpm as double)                     as tpm,
    cast(log2_tpm_plus1 as double)          as log2_tpm_plus1,
    qc_pass::boolean                        as qc_pass,
    nullif(trim(qc_fail_reason), '')        as qc_fail_reason,
    current_timestamp()                     as dbt_updated_at

from source
where
    gene_id    is not null
    and sample_id  is not null
    and tpm        >= 0
    and raw_count  >= 0
