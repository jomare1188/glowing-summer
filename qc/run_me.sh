#!/bin/bash
# conda activate qc_illumina   (fastp, fastqc, multiqc, seqkit, kraken2, jq, parallel)

KRAKEN_DB=/dados04/jorge/databases/kraken2/k2_pluspfp_16_GB_20260626 \
    bash qc_pipeline_v2.sh ../raw/ qc_results/ 200
