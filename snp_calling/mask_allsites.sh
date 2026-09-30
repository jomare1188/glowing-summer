#!/bin/bash
################################################################################
# MASKED ALL-SITES VCFs  (feedback.txt item 2; base for individual
# heterozygosity and ROH)
#
# Input: the unmasked all-sites VCFs of generate_all_vcf.sh (variant +
# invariant; SNP hard filters; --max-missing 0.8 and meanDP bounds on both;
# no MAF/HWE/LD). They are left untouched. Per chromosome, in parallel:
#
#   1. exclude mask_paralog_union.bed (v2)
#   2. genotype filters, identical for variant and invariant calls:
#        FORMAT/DP < GT_MIN_DP (3)                         -> ./.
#        het call failing binom(FORMAT/AD) < GT_BINOM_P (0.01) -> ./.
#      (FORMAT/DP, not GQ: invariant sites carry RGQ, not GQ)
#   3. AC/AN recomputed; sites with > 20 % missing genotypes dropped again
#   Sites left monomorphic by step 2 stay in the file: they are still
#   callable sites for the heterozygosity denominator.
#
# Why these genotype filters: tuning table in results/het_checks/
# gt_filter_tuning.tsv and in feedback.txt / README. In short: at 20x only
# 2-4 % of true het calls should look unbalanced, but 33-58 % did; the
# binomial test at p < 0.01 removes most of them and is not over-strict;
# DP < 3 (not 5) keeps 4 % more SNPs in d20 and 40 % more in d09 at the same
# quality.
#
# Outputs, per dataset, in results/<tag>/allsites/:
#   cohort.allsites.masked_v2.vcf.gz (+ .tbi)
#   allsites_masked_summary.tsv     per chromosome: invariant, SNP records,
#                                   polymorphic SNPs, sites total
#   per_tree_heterozygosity.tsv     per tree: callable sites, het sites, H
#
# Usage:  bash mask_allsites.sh [d20|d09|all]           (SNP_call env active)
#         TEST_REGION=Chr17:1-2000000 bash mask_allsites.sh d20
################################################################################

set -euo pipefail

WORK_DIR="/dados04/jorge/CALLICARPA_2/snp_calling"
REF_FAI="/dados04/jorge/CALLICARPA_2/CallicarpaGenome/car_asm.fa.fai"
MASK="${WORK_DIR}/results/masks/mask_paralog_union.bed"
GT_MIN_DP=${GT_MIN_DP:-3}
GT_BINOM_P=${GT_BINOM_P:-0.01}
MAX_MISSING_FRAC=0.2            # = --max-missing 0.8
JOBS=${JOBS:-18}
TEST_REGION="${TEST_REGION:-}"
WHICH="${1:-all}"

log() { echo "[$(date +'%F %T')] $*"; }
command -v bcftools >/dev/null || { echo "conda activate SNP_call" >&2; exit 1; }

# One region: mask -> genotype filters -> missingness -> piece + per-tree counts
piece() {
    set -euo pipefail
    local in=$1 name=$2 reg=$3 pd=$4
    local out="${pd}/${name}.vcf.gz"
    local rflag=(-r "${reg}"); [[ -f "${reg}" ]] && rflag=(-R "${reg}")
    bcftools view "${rflag[@]}" -T "^${MASK}" -Ou "${in}" \
      | bcftools +setGT -Ou -- -t q -n . -i "FMT/DP<${GT_MIN_DP}" 2> /dev/null \
      | bcftools +setGT -Ou -- -t q -n . -i "GT=\"het\" && binom(FMT/AD)<${GT_BINOM_P}" 2> /dev/null \
      | bcftools +fill-tags -Ou -- -t AC,AN 2> /dev/null \
      | bcftools view -i "F_MISSING<=${MAX_MISSING_FRAC}" -Oz -o "${out}"
    bcftools index -t -f "${out}"
    # counts: invariant / SNP records / polymorphic; per tree callable and het
    bcftools query -f '%ALT\t%INFO/AC\t%INFO/AN[\t%GT]\n' "${out}" \
      | awk -F'\t' -v OFS='\t' -v name="${name}" -v cnt="${pd}/${name}.counts" -v tree="${pd}/${name}.trees" '
          { if($1==".") inv++; else { snp++; n=split($2,ac,","); s=0; for(i=1;i<=n;i++) s+=ac[i]; if(s>0 && s<$3) poly++ }
            for(i=4;i<=NF;i++){ g=$i; if(g ~ /\./) continue; c[i]++; split(g,a,/[\/|]/); if(a[1]!=a[2]) h[i]++ }
            last=NF }
          END{ print name, inv+0, snp+0, poly+0, inv+snp > cnt
               for(i=4;i<=last;i++) print name, i-3, c[i]+0, h[i]+0 > tree }'
}
export -f piece
export MASK GT_MIN_DP GT_BINOM_P MAX_MISSING_FRAC

run_dataset() {
    local tag=$1
    local ad="${WORK_DIR}/results/${tag}/allsites"
    local in="${ad}/cohort.allsites.final.vcf.gz"
    local sfx=""; [[ -n "${TEST_REGION}" ]] && sfx=".test"
    local out="${ad}/cohort.allsites.masked_v2${sfx}.vcf.gz"
    local pd="${WORK_DIR}/tmp/mask_allsites_${tag}${sfx}"
    [[ -f "${in}.tbi" ]] || { echo "missing ${in}" >&2; exit 1; }
    rm -rf "${pd}"; mkdir -p "${pd}"
    log "${tag}: mask v2 + GT filters (DP < ${GT_MIN_DP}, het binom p < ${GT_BINOM_P}) + max-missing 0.8"

    # Regions in reference order: Chr01..Chr17, then all scaffolds as one job
    local jobs="${pd}/regions.tsv"
    if [[ -n "${TEST_REGION}" ]]; then
        printf 'test\t%s\n' "${TEST_REGION}" > "${jobs}"
    else
        awk '$1 ~ /^Chr/ {print $1"\t"$1}' "${REF_FAI}" > "${jobs}"
        awk -v OFS='\t' '$1 !~ /^Chr/ {print $1, 0, $2}' "${REF_FAI}" > "${pd}/scaffolds.bed"
        printf 'scaffolds\t%s\n' "${pd}/scaffolds.bed" >> "${jobs}"
    fi
    PARALLEL_SHELL=/bin/bash parallel --colsep '\t' --halt now,fail=1 -j "${JOBS}" \
        piece "${in}" {1} {2} "${pd}" :::: "${jobs}"

    local names; mapfile -t names < <(cut -f1 "${jobs}")
    local pieces=(); for n in "${names[@]}"; do pieces+=("${pd}/${n}.vcf.gz"); done
    bcftools concat --threads 8 -Oz -o "${out}" "${pieces[@]}"
    bcftools index -t -f "${out}"

    local summ="${ad}/allsites_masked_summary${sfx}.tsv"
    printf 'interval\tinvariant_sites\tsnp_records\tpolymorphic_snps\ttotal_sites\n' > "${summ}"
    for n in "${names[@]}"; do cat "${pd}/${n}.counts"; done >> "${summ}"
    awk -F'\t' -v OFS='\t' 'NR>1{for(i=2;i<=5;i++) t[i]+=$i} END{print "TOTAL", t[2], t[3], t[4], t[5]}' "${summ}" >> "${summ}"

    local het="${ad}/per_tree_heterozygosity${sfx}.tsv"
    bcftools query -l "${out}" > "${pd}/samples.txt"
    cat "${pd}"/*.trees | awk -F'\t' -v OFS='\t' 'NR==FNR{nm[FNR]=$1; next} {c[$2]+=$3; h[$2]+=$4; k=($2>k?$2:k)}
        END{ print "tree", "callable_sites", "het_sites", "H"
             for(i=1;i<=k;i++) printf "%s\t%d\t%d\t%.6f\n", nm[i], c[i], h[i], h[i]/c[i] }' "${pd}/samples.txt" - > "${het}"

    bcftools stats "${out}" > "${out%.vcf.gz}.stats.txt"
    rm -rf "${pd}"
    log "${tag}: $(tail -1 "${summ}" | cut -f2) invariant + $(tail -1 "${summ}" | cut -f3) SNP records ($(tail -1 "${summ}" | cut -f4) polymorphic) -> ${out}"
    column -t "${het}"
}

for t in d20 d09; do
    if [[ "${WHICH}" == all || "${WHICH}" == "${t}" ]]; then run_dataset "${t}"; fi
done
