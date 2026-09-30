#!/bin/bash
################################################################################
# MASKED SNP VCF FOR DELIVERY  (Tatiana, 30 Sep 2026)
#
# SNPs of one dataset with the paralog mask (v2) applied, and NO MAF, NO HWE
# and NO LD pruning. cohort.snps.final.vcf.gz cannot be used: it already went
# through --maf 0.01 --hwe 0.01. So this starts from the GATK hard-filtered
# PASS set and repeats the SAME missingness and depth filters as
# snp_datasets.sh (--max-missing 0.8, meanDP scaled to the target depth).
#
#   cohort.snps.pass.vcf.gz
#     -> exclude mask_paralog_union.bed, biallelic SNPs only
#     -> vcftools --remove-indels --max-missing 0.8 --min/max-meanDP
#        (--recode-INFO-all: INFO/DP, FORMAT/AD/DP/GQ/PL kept)
#     -> drop sites monomorphic in the cohort (AC = 0 or AC = AN). They are not
#        SNPs within these trees; with 9 trees this is the only thing a
#        MAF 0.01 filter would remove anyway (lowest possible MAF is 1/18).
#
# Usage: bash masked_snps.sh [d20|d09]     (SNP_call env active)
# Output: results/<tag>/masked/cohort.snps.masked_v2.vcf.gz (+ .tbi, stats)
################################################################################

set -euo pipefail

WORK_DIR="/dados04/jorge/CALLICARPA_2/snp_calling"
TAG="${1:-d20}"
case "${TAG}" in d09) TARGET=9 ;; d20) TARGET=20 ;; *) echo "use d09 or d20" >&2; exit 1 ;; esac
# Same bounds as scale_depth_filters() in snp_datasets.sh
MIN_MEAN_DP=$(awk -v t="${TARGET}" 'BEGIN{printf "%d", (t/2 < 3 ? 3 : t/2)}')
MAX_MEAN_DP=$(awk -v t="${TARGET}" 'BEGIN{printf "%d", t*3}')
VCFTOOLS_BIN="${VCFTOOLS_BIN:-/home/genomics/miniconda3/envs/popgen_tools/bin/vcftools}"

PASS="${WORK_DIR}/results/${TAG}/filtered/cohort.snps.pass.vcf.gz"
MASK="${WORK_DIR}/results/masks/mask_paralog_union.bed"
OUTD="${WORK_DIR}/results/${TAG}/masked"
OUT="${OUTD}/cohort.snps.masked_v2.vcf.gz"
tmp="${WORK_DIR}/tmp/masked_snps_${TAG}"

command -v bcftools >/dev/null || { echo "conda activate SNP_call" >&2; exit 1; }
mkdir -p "${OUTD}" "${tmp}"
echo "[$(date +'%F %T')] ${TAG}: mask $(basename "${MASK}"), max-missing 0.8, meanDP ${MIN_MEAN_DP}-${MAX_MEAN_DP}, no MAF/HWE/LD"

bcftools view -T "^${MASK}" -m2 -M2 -v snps -Oz -o "${tmp}/masked.vcf.gz" "${PASS}"
"${VCFTOOLS_BIN}" --gzvcf "${tmp}/masked.vcf.gz" \
    --remove-indels --max-missing 0.8 \
    --min-meanDP "${MIN_MEAN_DP}" --max-meanDP "${MAX_MEAN_DP}" \
    --recode --recode-INFO-all --stdout 2> "${OUTD}/vcftools.log" \
  | bcftools view -e 'INFO/AC==0 || INFO/AC==INFO/AN' -Oz -o "${OUT}"
tabix -f -p vcf "${OUT}"
bcftools stats "${OUT}" > "${OUTD}/cohort.snps.masked_v2.stats.txt"
rm -rf "${tmp}"
echo "[$(date +'%F %T')] ${TAG}: $(bcftools index -n "${OUT}") SNPs -> ${OUT}"
