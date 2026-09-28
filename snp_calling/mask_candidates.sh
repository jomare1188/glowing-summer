#!/bin/bash
################################################################################
# STRICTER MASK CANDIDATES  (feedback.txt, Tatiana 28 Sep 2026, check 3)
#
# het_checks.sh showed that paralogs remain after mask_paralog_union.bed: on
# masked d20, SNPs with all 9 trees heterozygous are 55x the HWE expectation,
# and removing the three high-Ho trees does not help. Each candidate below is
# the current union PLUS one extra layer, scored with het_checks.sh on d20.
# All inputs already exist; nothing is re-called.
#
#   ab      pooled alt-read fraction across het calls outside 0.3-0.7
#           (native13 full depth, >= 20 het reads), GATK "ABHet"-style
#   eh13    ExcessHet >= 13 (p < 0.05 at n = 13) instead of >= 20
#   h75     >= 75% of called trees heterozygous (native13, >= 10 trees called)
#   win1kb  paralog sites (ExcessHet >= 20) extended by 1 kb each side
#   win500  the same, 500 bp each side
#   win1kb_h75  win1kb plus h75
#
# Caution: eh13 and h75 are built from heterozygosity itself, so they push Fis
# toward 0 by construction. Judge them on the spectrum and SNPs lost too.
#
# Usage: bash mask_candidates.sh     (SNP_call env active)
# Output: results/masks/candidates/*.bed, results/het_checks/mask_candidates.tsv
################################################################################

set -euo pipefail

WORK_DIR="/dados04/jorge/CALLICARPA_2/snp_calling"
MASKS="${WORK_DIR}/results/masks"
CAND="${MASKS}/candidates"
STATS="${MASKS}/paralog/paralog_site_stats.tsv.gz"
UNION="${MASKS}/mask_paralog_union.bed"
D20="${WORK_DIR}/results/d20/filtered/cohort.snps.final.vcf.gz"
BT="/home/genomics/miniconda3/envs/paralog_masks/bin/bedtools"
G="${MASKS}/genome.txt"
HC="${WORK_DIR}/results/het_checks"

command -v bcftools >/dev/null || { echo "conda activate SNP_call" >&2; exit 1; }
mkdir -p "${CAND}"

# Extra layer from per-site stats: $1 name, $2 awk condition on the stats columns
# (3 excesshet, 4 n_called, 6 H, 7 het_reads, 8 D).
layer_from_stats() {
    zcat "${STATS}" | awk -F'\t' -v OFS='\t' 'NR>1 && ('"$2"') {print $1,$2-1,$2}' \
      | "${BT}" sort -g "${G}" -i - | "${BT}" merge -i - > "${CAND}/extra_$1.bed"
}
# D = (ref - n/2)/sqrt(n/4)  =>  ref fraction = 0.5 + D/(2 sqrt(n)); alt = 1 - ref
layer_from_stats ab   '$7>=20 && $8!="NA" && (0.5+$8/(2*sqrt($7)) < 0.3 || 0.5+$8/(2*sqrt($7)) > 0.7)'
layer_from_stats eh13 '$3>=13'
layer_from_stats h75  '$4>=10 && $6>=0.75'
"${BT}" slop -b 1000 -g "${G}" -i "${MASKS}/mask_paralog_sites.bed" \
  | "${BT}" sort -g "${G}" -i - | "${BT}" merge -i - > "${CAND}/extra_win1kb.bed"
"${BT}" slop -b 500 -g "${G}" -i "${MASKS}/mask_paralog_sites.bed" \
  | "${BT}" sort -g "${G}" -i - | "${BT}" merge -i - > "${CAND}/extra_win500.bed"
cat "${CAND}/extra_win1kb.bed" "${CAND}/extra_h75.bed" \
  | "${BT}" sort -g "${G}" -i - | "${BT}" merge -i - > "${CAND}/extra_win1kb_h75.bed"

tab="${HC}/mask_candidates.tsv"
printf 'candidate\textra_bp\tunion_bp\tpct_genome\td20_snps_kept\tall9het_obs\tall9het_exp\tall9het_obs_over_exp\t8of9_obs_over_exp\tFis_site\tFis_corrected\n' > "${tab}"

score() {  # $1 candidate name, $2 union bed
    local c=$1 u=$2 tag="d20_mask_$1"
    VCF="${D20}" MASK="${u}" TAG="${tag}" bash "${WORK_DIR}/het_checks.sh" > "${WORK_DIR}/logs/het_${tag}.log" 2>&1
    local ubp; ubp=$(awk '{s+=$3-$2} END{print s}' "${u}")
    local xbp=0; [[ -f "${CAND}/extra_${c}.bed" ]] && xbp=$(awk '{s+=$3-$2} END{print s}' "${CAND}/extra_${c}.bed")
    local sp="${HC}/het_spectrum_${tag}.tsv"
    local f; f=$(awk -F'\t' -v t="${tag}" '$1==t' "${HC}/fis_summary.tsv")
    printf '%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${c}" "${xbp}" "${ubp}" \
        "$(awk -v b="${ubp}" 'BEGIN{print 100*b/506362408}')" "$(cut -f3 <<< "${f}")" \
        "$(awk -F'\t' '$1==9{print $2}' "${sp}")" "$(awk -F'\t' '$1==9{print $3}' "${sp}")" \
        "$(awk -F'\t' '$1==9{print $4}' "${sp}")" "$(awk -F'\t' '$1==8{print $4}' "${sp}")" \
        "$(cut -f5 <<< "${f}")" "$(cut -f8 <<< "${f}")" >> "${tab}"
    echo "  scored ${c}"
}

score current "${UNION}"
for c in ab eh13 h75 win1kb win500 win1kb_h75; do
    cat "${UNION}" "${CAND}/extra_${c}.bed" | "${BT}" sort -g "${G}" -i - | "${BT}" merge -i - > "${CAND}/union_plus_${c}.bed"
    score "${c}" "${CAND}/union_plus_${c}.bed"
done
column -t "${tab}"
