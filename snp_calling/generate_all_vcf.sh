#!/bin/bash
################################################################################
# ALL-SITES VCFs (variant + invariant) FOR THE DEPTH-EQUALISED DATASETS
#
# Companion to snp_datasets.sh, which is left untouched. That pipeline emits
# SNPs only; diversity (pi) and per-individual heterozygosity also need the
# invariant sites, otherwise the denominator is lost and estimates inflate.
#
# Input: the per-sample GVCFs retained by snp_datasets.sh
#        (results/{d09,d20}/gvcf/*.g.vcf.gz, -ERC GVCF, reference blocks
#        included), so nothing is re-aligned or re-called.
#
# Per dataset, per interval (Chr01..Chr17 + one job for all scaffolds), in
# parallel:
#   CombineGVCFs -> GenotypeGVCFs --include-non-variant-sites
#   split: invariant (ALT=".")  |  SNPs -> GATK hard filters -> PASS
#   SAME population filters on both parts, taken from snp_datasets.sh:
#       --remove-indels --max-missing 0.8 --min-meanDP/--max-meanDP (scaled)
#   NO --maf, NO --hwe, NO LD pruning
#   merge back into one position-sorted VCF
# then the per-interval pieces are concatenated in reference order.
#
# Chromosome names stay as in the reference (Chr01...): PLINK is never used.
#
# Usage:  bash generate_all_vcf.sh [d09|d20|all]            (default: all)
#         TEST_INTERVAL=Chr17:1-2000000 bash generate_all_vcf.sh d20
#             -> runs one region only, into results/<tag>/allsites_test/
################################################################################

set -euo pipefail

################################################################################
# CONFIGURATION  (filters copied verbatim from snp_datasets.sh)
################################################################################

PROJECT_DIR="/dados04/jorge/CALLICARPA_2"
REF_GENOME="${PROJECT_DIR}/CallicarpaGenome/car_asm.fa"
WORK_DIR="${PROJECT_DIR}/snp_calling"

OUT_DIR="${WORK_DIR}/results"
LOG_DIR="${WORK_DIR}/logs"
TMP_DIR="${WORK_DIR}/tmp"
CHECKPOINT_DIR="${WORK_DIR}/checkpoints"
MAIN_LOG="${LOG_DIR}/allsites_main.log"

# Datasets: <tag>:<target depth in X>
DATASETS=("d09:9" "d20:20")

# Resources. One job = one interval (CombineGVCFs + GenotypeGVCFs + filters).
MAX_PARALLEL_JOBS=${MAX_PARALLEL_JOBS:-18}
MEM_JOB="${MEM_JOB:-16G}"

# GATK hard filters (site level; SNPs only -- invariant sites carry no QD/FS/MQ)
HARD_FILTERS=(
    "QD_filter:QD < 2.0"
    "FS_filter:FS > 60.0"
    "MQ_filter:MQ < 40.0"
    "SOR_filter:SOR > 3.0"
    "MQRankSum_filter:MQRankSum < -12.5"
    "ReadPosRankSum_filter:ReadPosRankSum < -8.0"
)

# Population filters applied identically to variant and invariant sites.
# Deliberately no MAF, no HWE: they remove real rare variants and real
# heterozygotes, which is exactly what diversity/heterozygosity measure.
POP_MAX_MISSING=0.8

DISK_WARN_PCT=85
DISK_ABORT_PCT=96

FORCE_RERUN=${FORCE_RERUN:-false}
TEST_INTERVAL="${TEST_INTERVAL:-}"

WHICH="${1:-all}"

################################################################################
# LOGGING AND CHECKPOINTS  (same helpers as snp_datasets.sh)
################################################################################

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

_emit() { echo -e "$1"; echo -e "$1" | sed 's/\x1b\[[0-9;]*m//g' >> "${MAIN_LOG}"; }
log()   { _emit "${GREEN}[$(date +'%F %T')]${NC} $1"; }
info()  { _emit "${BLUE}[INFO]${NC} $1"; }
warn()  { _emit "${YELLOW}[WARN]${NC} $1"; }
error() { _emit "${RED}[ERROR]${NC} $1"; exit 1; }

checkpoint_exists() { [[ -f "${CHECKPOINT_DIR}/$1.$2.done" ]] && ! $FORCE_RERUN; }
create_checkpoint() { date +%s > "${CHECKPOINT_DIR}/$1.$2.done"; }

disk_pct() { df -P "${WORK_DIR}" | awk 'NR==2 {gsub(/%/,"",$5); print $5}'; }

check_disk() {
    local pct; pct=$(disk_pct)
    (( pct >= DISK_ABORT_PCT )) && error "disk ${pct}% full, aborting to avoid a truncated write"
    (( pct >= DISK_WARN_PCT )) && warn "disk ${pct}% full"
    return 0
}

# Same bounds as snp_datasets.sh: half and three times the target depth.
scale_depth_filters() {
    local target=$1
    MIN_MEAN_DP=$(awk -v t="${target}" 'BEGIN{printf "%d", (t/2 < 3 ? 3 : t/2)}')
    MAX_MEAN_DP=$(awk -v t="${target}" 'BEGIN{printf "%d", t*3}')
}

################################################################################
# PREFLIGHT
################################################################################

preflight() {
    mkdir -p "${LOG_DIR}" "${TMP_DIR}" "${CHECKPOINT_DIR}"
    touch "${MAIN_LOG}"
    log "All-sites VCFs — preflight"

    local missing=()
    for t in gatk bcftools tabix bgzip parallel vcftools; do
        command -v "$t" >/dev/null || missing+=("$t")
    done
    if (( ${#missing[@]} > 0 )); then
        error "missing tools: ${missing[*]}
  Activate both environments:
    source /home/genomics/miniconda3/etc/profile.d/conda.sh
    conda activate SNP_call
    export PATH=\"/home/genomics/miniconda3/envs/popgen_tools/bin:\$PATH\""
    fi
    [[ -f "${REF_GENOME}.fai" ]] || error "reference index not found: ${REF_GENOME}.fai"
    [[ -n "${TEST_INTERVAL}" ]] && warn "TEST MODE: only ${TEST_INTERVAL}"
    return 0
}

################################################################################
# PER-INTERVAL WORK  (run under GNU parallel)
################################################################################

# $1 tag  $2 interval name (file-safe)  $3 interval (-L value)  $4 allsites dir
# $5 ckpt prefix  $6 minDP  $7 maxDP
process_interval() {
    set -euo pipefail
    local tag=$1 name=$2 iv=$3 ad=$4 ck=$5 min_dp=$6 max_dp=$7
    local pd="${ad}/pieces" tmp="${TMP_DIR}/${ck}_${tag}_${name}"
    local combined="${pd}/${name}.g.vcf.gz"
    local raw="${pd}/${name}.raw.vcf.gz"
    local snps="${pd}/${name}.snps.vcf.gz"
    local marked="${pd}/${name}.snps.marked.vcf.gz"
    local pass="${pd}/${name}.snps.pass.vcf.gz"
    local inv="${pd}/${name}.invariant.vcf.gz"
    local snps_f="${pd}/${name}.snps.filt.vcf.gz"
    local inv_f="${pd}/${name}.invariant.filt.vcf.gz"
    local final="${pd}/${name}.allsites.filt.vcf.gz"
    local lg="${LOG_DIR}/${ck}_${tag}_${name}"

    if checkpoint_exists "${ck}_${tag}" "${name}"; then
        info "  ${tag}/${name}: done, skipping"
        return 0
    fi
    check_disk
    mkdir -p "${tmp}"

    if [[ ! -f "${raw}.tbi" ]]; then
        log "  ${tag}/${name}: CombineGVCFs"
        local args=()
        while read -r s; do args+=(-V "${OUT_DIR}/${tag}/gvcf/${s}.g.vcf.gz"); done \
            < "${OUT_DIR}/metrics/samples_${tag}.txt"
        gatk --java-options "-Xmx${MEM_JOB}" CombineGVCFs \
            -R "${REF_GENOME}" "${args[@]}" -L "${iv}" -O "${combined}" \
            --tmp-dir "${tmp}" &> "${lg}_combine.log"

        log "  ${tag}/${name}: GenotypeGVCFs --include-non-variant-sites"
        gatk --java-options "-Xmx${MEM_JOB}" GenotypeGVCFs \
            -R "${REF_GENOME}" -V "${combined}" -L "${iv}" -O "${raw}" \
            --include-non-variant-sites \
            --tmp-dir "${tmp}" &> "${lg}_genotype.log"
        rm -f "${combined}" "${combined}.tbi"
    fi

    # Invariant sites: no ALT allele at all.
    bcftools view --max-alleles 1 -O z -o "${inv}" "${raw}"
    tabix -f -p vcf "${inv}"

    # Variant sites: same selection and hard filters as snp_datasets.sh.
    gatk SelectVariants -R "${REF_GENOME}" -V "${raw}" \
        -select-type SNP -O "${snps}" &> "${lg}_select.log"
    local fargs=()
    for f in "${HARD_FILTERS[@]}"; do fargs+=(--filter-name "${f%%:*}" --filter-expression "${f#*:}"); done
    gatk VariantFiltration -R "${REF_GENOME}" -V "${snps}" -O "${marked}" \
        "${fargs[@]}" &> "${lg}_filtration.log"
    bcftools view -f PASS -O z -o "${pass}" "${marked}"
    tabix -f -p vcf "${pass}"

    # Identical population filters on both parts; no --maf, no --hwe.
    local part in out
    for part in snps invariant; do
        if [[ ${part} == snps ]]; then in="${pass}"; out="${snps_f}"; else in="${inv}"; out="${inv_f}"; fi
        vcftools --gzvcf "${in}" \
            --remove-indels \
            --max-missing "${POP_MAX_MISSING}" \
            --min-meanDP "${min_dp}" \
            --max-meanDP "${max_dp}" \
            --recode --recode-INFO-all --stdout 2> "${lg}_popfilter_${part}.log" \
          | bgzip > "${out}"
        tabix -f -p vcf "${out}"
    done

    bcftools concat -a -O u "${snps_f}" "${inv_f}" \
      | bcftools sort -T "${tmp}" -O z -o "${final}" 2> "${lg}_sort.log"
    tabix -f -p vcf "${final}"

    # Per-interval counts: raw, invariant kept, SNPs kept
    printf '%s\t%s\t%s\t%s\n' "${name}" \
        "$(bcftools view -H "${raw}" | wc -l)" \
        "$(bcftools index -n "${inv_f}")" \
        "$(bcftools index -n "${snps_f}")" > "${pd}/${name}.counts.tsv"

    rm -f "${snps}" "${snps}.tbi" "${marked}" "${marked}.tbi" "${pass}" "${pass}.tbi" \
          "${inv}" "${inv}.tbi" "${snps_f}" "${snps_f}.tbi" "${inv_f}" "${inv_f}.tbi"
    rm -rf "${tmp}"
    create_checkpoint "${ck}_${tag}" "${name}"
    log "  ${tag}/${name}: done"
}
export -f process_interval _emit log info warn error checkpoint_exists create_checkpoint disk_pct check_disk
export OUT_DIR LOG_DIR TMP_DIR CHECKPOINT_DIR MAIN_LOG REF_GENOME WORK_DIR MEM_JOB \
       POP_MAX_MISSING FORCE_RERUN DISK_WARN_PCT DISK_ABORT_PCT RED GREEN YELLOW BLUE NC
# Arrays cannot be exported; hand HARD_FILTERS to the workers through a
# declare string that each worker evaluates.
export HARD_FILTERS_DECL
HARD_FILTERS_DECL=$(declare -p HARD_FILTERS)

################################################################################
# PER DATASET
################################################################################

build_dataset() {
    local tag=$1 target=$2
    local ad ck
    if [[ -n "${TEST_INTERVAL}" ]]; then
        ad="${OUT_DIR}/${tag}/allsites_test"; ck="allsites_test"
    else
        ad="${OUT_DIR}/${tag}/allsites"; ck="allsites"
    fi
    local pd="${ad}/pieces"
    local final="${ad}/cohort.allsites.final.vcf.gz"
    local raw_all="${ad}/raw/cohort.allsites.raw.vcf.gz"
    mkdir -p "${pd}" "${ad}/raw"

    scale_depth_filters "${target}"
    log "Dataset ${tag} (${target}x, $(wc -l < "${OUT_DIR}/metrics/samples_${tag}.txt") samples): meanDP ${MIN_MEAN_DP}-${MAX_MEAN_DP}, max-missing ${POP_MAX_MISSING}, no MAF/HWE/LD"

    if checkpoint_exists "${ck}_${tag}" merged; then
        info "  ${tag}: already merged, skipping"
        return 0
    fi

    # Jobs file: <name>\t<-L value>, in reference order.
    local jobs="${ad}/intervals.tsv"
    if [[ -n "${TEST_INTERVAL}" ]]; then
        printf '%s\t%s\n' "$(tr ':-' '__' <<< "${TEST_INTERVAL}")" "${TEST_INTERVAL}" > "${jobs}"
    else
        awk '$1 ~ /^Chr/ {print $1}' "${REF_GENOME}.fai" | awk '{print $1"\t"$1}' > "${jobs}"
        awk '$1 !~ /^Chr/ {print $1}' "${REF_GENOME}.fai" > "${ad}/scaffolds.list"
        printf 'scaffolds\t%s\n' "${ad}/scaffolds.list" >> "${jobs}"
    fi
    info "  ${tag}: $(wc -l < "${jobs}") intervals, ${MAX_PARALLEL_JOBS} in parallel"

    # shellcheck disable=SC2016
    PARALLEL_SHELL=/bin/bash parallel --colsep '\t' -j "${MAX_PARALLEL_JOBS}" --halt now,fail=1 \
        'eval "${HARD_FILTERS_DECL}"; process_interval '"${tag}"' {1} {2} '"${ad} ${ck} ${MIN_MEAN_DP} ${MAX_MEAN_DP}" \
        :::: "${jobs}"

    log "  ${tag}: concatenating $(wc -l < "${jobs}") pieces in reference order"
    local names; mapfile -t names < <(cut -f1 "${jobs}")
    local fin=() rawp=()
    for n in "${names[@]}"; do fin+=("${pd}/${n}.allsites.filt.vcf.gz"); rawp+=("${pd}/${n}.raw.vcf.gz"); done
    bcftools concat --threads 8 -O z -o "${final}" "${fin[@]}"
    tabix -f -p vcf "${final}"
    bcftools concat --threads 8 -O z -o "${raw_all}" "${rawp[@]}"
    tabix -f -p vcf "${raw_all}"

    # Summary table
    local summ="${ad}/allsites_summary.tsv"
    printf 'interval\traw_sites\tinvariant_kept\tsnps_kept\n' > "${summ}"
    for n in "${names[@]}"; do cat "${pd}/${n}.counts.tsv"; done >> "${summ}"
    awk -F'\t' 'NR>1{r+=$2;i+=$3;s+=$4} END{printf "TOTAL\t%d\t%d\t%d\n",r,i,s}' "${summ}" >> "${summ}"
    bcftools stats "${final}" > "${ad}/cohort.allsites.final.stats.txt"

    for f in "${fin[@]}" "${rawp[@]}"; do rm -f "${f}" "${f}.tbi"; done
    create_checkpoint "${ck}_${tag}" merged

    local tot; tot=$(tail -1 "${summ}")
    log "  ${tag}: $(cut -f3 <<< "${tot}") invariant + $(cut -f4 <<< "${tot}") SNPs kept of $(cut -f2 <<< "${tot}") raw sites -> ${final}"
}

################################################################################
# MAIN
################################################################################

main() {
    preflight
    local ran=0
    for tag_target in "${DATASETS[@]}"; do
        local tag="${tag_target%%:*}" target="${tag_target##*:}"
        [[ "${WHICH}" == all || "${WHICH}" == "${tag}" ]] || continue
        build_dataset "${tag}" "${target}"
        ran=1
    done
    (( ran )) || error "unknown dataset '${WHICH}' (use: d09|d20|all)"
    info "  disk $(disk_pct)% used, $(df -Ph "${WORK_DIR}" | awk 'NR==2{print $4}') free"
    log "All-sites VCFs finished"
}

main
