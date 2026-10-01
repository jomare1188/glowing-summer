#!/bin/bash
################################################################################
# MASKED SNP VCF  (feedback.txt item 3; consistent with items 1 and 2)
#
# The SNP VCF is taken straight from the masked all-sites VCF built by
# mask_allsites.sh, so both have exactly the same sites and genotypes:
#   paralog mask v2; GATK hard filters; --max-missing 0.8 and meanDP bounds;
#   genotype filters DP < 3 and het binomial balance p < 0.01; max-missing 0.8
#   re-applied. No MAF, no HWE, no LD pruning.
# Kept: biallelic SNPs polymorphic in the cohort (0 < AC < AN). Dropped:
# invariant sites, multiallelic/spanning-deletion records, and SNP records
# left monomorphic (fixed ALT, or no ALT carrier after the genotype filters).
#
# History: the 30 Sep 2026 version started from the PASS set without the
# genotype filters (d20 401,017 SNPs). It was not consistent with item 2 and
# is kept in results/<tag>/masked/old_provisional_20260930/.
#
# Usage: bash masked_snps.sh [d20|d09]     (SNP_call env active)
# Output: results/<tag>/masked/cohort.snps.masked_v2_gtfilt.vcf.gz (+ .tbi, stats)
################################################################################

set -euo pipefail

WORK_DIR="/dados04/jorge/CALLICARPA_2/snp_calling"
TAG="${1:-d20}"
case "${TAG}" in d09|d20) ;; *) echo "use d09 or d20" >&2; exit 1 ;; esac

IN="${WORK_DIR}/results/${TAG}/allsites/cohort.allsites.masked_v2.vcf.gz"
OUTD="${WORK_DIR}/results/${TAG}/masked"
OUT="${OUTD}/cohort.snps.masked_v2_gtfilt.vcf.gz"

command -v bcftools >/dev/null || { echo "conda activate SNP_call" >&2; exit 1; }
[[ -f "${IN}.tbi" ]] || { echo "missing ${IN}; run mask_allsites.sh ${TAG}" >&2; exit 1; }
mkdir -p "${OUTD}"
echo "[$(date +'%F %T')] ${TAG}: polymorphic biallelic SNPs from $(basename "${IN}")"

bcftools view --threads 4 -m2 -M2 -v snps -i 'INFO/AC>0 && INFO/AC<INFO/AN' -Oz -o "${OUT}" "${IN}"
tabix -f -p vcf "${OUT}"
bcftools stats "${OUT}" > "${OUT%.vcf.gz}.stats.txt"
echo "[$(date +'%F %T')] ${TAG}: $(bcftools index -n "${OUT}") SNPs -> ${OUT}"
