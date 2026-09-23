#!/bin/bash
################################################################################
# NATIVE-DEPTH BAMs AND JOINT CALLS  (feedback.txt step 1; also feeds steps 4, 6)
#
# snp_datasets.sh deletes its BAMs and only ever calls DOWNSAMPLED reads. The
# paralog masks need data at full depth:
#   - max-depth mask      -> native BAMs (mosdepth, in paralog_masks.sh)
#   - ExcessHet / H / D   -> joint calls made from the native BAMs
# and step 4 (clonality, ANGSD/NgsRelate) needs the native BAMs of all 18
# samples. snp_datasets.sh is left untouched; alignment commands are copied
# verbatim from its process_sample().
#
# Phases:
#   align  all 18 samples: bwa mem -> sort -> MarkDuplicates; markdup BAM KEPT
#          (+ .bai), flagstat, duplicate metrics, usable depth (-F 0xD04)
#   call   HaplotypeCaller -ERC GVCF on the 13 d09 samples, scattered per
#          sample x interval (Chr01..Chr17 + one job for all scaffolds)
#   joint  per interval CombineGVCFs -> GenotypeGVCFs (variant sites), then
#          concat -> results/native/raw/native13.raw.vcf.gz
#
# Why 13 and not 18 for calling: ExcessHet and heterozygosity statistics assume
# independent individuals with adequate depth. 563/564 are known duplicates of
# 562/561, and 553/557/558 sit at 3-4x where heterozygotes drop out. The 13 are
# exactly the d09 members (>= 9x usable depth).
#
# Usage:  bash native_depth.sh [all|align|call|joint]        (default: all)
#         TEST_INTERVAL=Chr17:1-2000000 bash native_depth.sh call   (+ joint)
#             -> one region only, into results/native_test/
################################################################################

set -euo pipefail

################################################################################
# CONFIGURATION
################################################################################

PROJECT_DIR="/dados04/jorge/CALLICARPA_2"
READS_DIR="${PROJECT_DIR}/qc/qc_results/trimmed_reads"
REF_GENOME="${PROJECT_DIR}/CallicarpaGenome/car_asm.fa"
WORK_DIR="${PROJECT_DIR}/snp_calling"
GENOME_SIZE=506362408

LOG_DIR="${WORK_DIR}/logs"
TMP_DIR="${WORK_DIR}/tmp"
CHECKPOINT_DIR="${WORK_DIR}/checkpoints"
MAIN_LOG="${LOG_DIR}/native_main.log"

TEST_INTERVAL="${TEST_INTERVAL:-}"
if [[ -n "${TEST_INTERVAL}" ]]; then
    NAT_DIR="${WORK_DIR}/results/native_test"; CK="native_test"
else
    NAT_DIR="${WORK_DIR}/results/native"; CK="native"
fi
BAM_DIR="${WORK_DIR}/results/native/bam"          # BAMs are shared with test mode
MET_DIR="${WORK_DIR}/results/native/metrics"
V2_DEPTH="${WORK_DIR}/results/metrics/depth_summary.tsv"
CALL_SAMPLES="${CALL_SAMPLES:-${WORK_DIR}/results/metrics/samples_d09.txt}"

# Resources
MAX_PARALLEL_SAMPLES=${MAX_PARALLEL_SAMPLES:-6}    # alignment jobs
THREADS_PER_SAMPLE=${THREADS_PER_SAMPLE:-10}       # bwa threads (as in snp_datasets.sh)
MEM_PER_SAMPLE="${MEM_PER_SAMPLE:-40G}"            # MarkDuplicates
HC_JOBS=${HC_JOBS:-40}                             # sample x interval jobs
HC_THREADS=${HC_THREADS:-4}
MEM_HC="${MEM_HC:-8G}"
JOINT_JOBS=${JOINT_JOBS:-18}
MEM_JOINT="${MEM_JOINT:-16G}"

MIN_BQ=20                  # as in snp_datasets.sh

DISK_WARN_PCT=85
DISK_ABORT_PCT=96
FORCE_RERUN=${FORCE_RERUN:-false}

PHASE="${1:-all}"

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

################################################################################
# PREFLIGHT
################################################################################

preflight() {
    mkdir -p "${LOG_DIR}" "${TMP_DIR}" "${CHECKPOINT_DIR}" "${BAM_DIR}" "${MET_DIR}" \
             "${NAT_DIR}"/{gvcf,pieces,raw}
    touch "${MAIN_LOG}"
    log "Native depth — preflight (phase ${PHASE})"

    local missing=()
    for t in bwa samtools gatk bcftools tabix bgzip parallel; do
        command -v "$t" >/dev/null || missing+=("$t")
    done
    (( ${#missing[@]} == 0 )) || error "missing tools: ${missing[*]} (conda activate SNP_call)"
    [[ -f "${REF_GENOME}.bwt" && -f "${REF_GENOME}.fai" ]] || error "reference indexes missing"
    [[ -f "${CALL_SAMPLES}" ]] || error "sample list not found: ${CALL_SAMPLES}"

    SAMPLES=()
    shopt -s nullglob
    for r1 in "${READS_DIR}"/*_R1_trimmed.fastq.gz; do
        local s; s=$(basename "${r1}" _R1_trimmed.fastq.gz)
        [[ -f "${READS_DIR}/${s}_R2_trimmed.fastq.gz" ]] || error "missing mate for ${r1}"
        SAMPLES+=("${s}")
    done
    shopt -u nullglob
    (( ${#SAMPLES[@]} > 0 )) || error "no trimmed reads in ${READS_DIR}"
    info "  ${#SAMPLES[@]} samples with reads; calling on $(wc -l < "${CALL_SAMPLES}") ($(basename "${CALL_SAMPLES}"))"
    [[ -n "${TEST_INTERVAL}" ]] && warn "TEST MODE: call/joint restricted to ${TEST_INTERVAL}"
    log "  disk $(disk_pct)% used, $(df -Ph "${WORK_DIR}" | awk 'NR==2{print $4}') free"
    return 0
}

################################################################################
# PHASE 1 — ALIGN (commands copied from snp_datasets.sh process_sample)
################################################################################

align_sample() {
    set -euo pipefail
    local sample=$1
    local r1="${READS_DIR}/${sample}_R1_trimmed.fastq.gz"
    local r2="${READS_DIR}/${sample}_R2_trimmed.fastq.gz"
    local sorted="${BAM_DIR}/${sample}.sorted.bam"
    local markdup="${BAM_DIR}/${sample}.markdup.bam"

    if checkpoint_exists "native_${sample}" alignment; then
        info "  ${sample}: alignment done, skipping"
    else
        log "  ${sample}: aligning"
        check_disk
        bwa mem -t "${THREADS_PER_SAMPLE}" -M \
            -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}\tPU:unit1" \
            "${REF_GENOME}" "${r1}" "${r2}" 2> "${LOG_DIR}/native_${sample}_bwa.log" \
          | samtools sort -@ 4 -m 2G -l 9 -T "${TMP_DIR}/native_${sample}.sort" -o "${sorted}" - \
            2> "${LOG_DIR}/native_${sample}_sort.log"
        samtools flagstat -@ 4 "${sorted}" > "${MET_DIR}/${sample}_flagstat.txt"
        create_checkpoint "native_${sample}" alignment
    fi

    if checkpoint_exists "native_${sample}" markdup; then
        info "  ${sample}: markdup done, skipping"
    else
        log "  ${sample}: marking duplicates"
        check_disk
        mkdir -p "${TMP_DIR}/native_${sample}"
        gatk --java-options "-Xmx${MEM_PER_SAMPLE}" MarkDuplicates \
            -I "${sorted}" -O "${markdup}" \
            -M "${MET_DIR}/${sample}_duplicate_metrics.txt" \
            --CREATE_INDEX false --VALIDATION_STRINGENCY SILENT \
            --COMPRESSION_LEVEL 9 --TMP_DIR "${TMP_DIR}/native_${sample}" \
            &> "${LOG_DIR}/native_${sample}_markdup.log"
        # samtools index gives the conventional <bam>.bai name (Picard writes <name>.bai)
        samtools index -@ 4 "${markdup}"
        rm -rf "${TMP_DIR}/native_${sample}"
        rm -f "${sorted}"
        create_checkpoint "native_${sample}" markdup
    fi

    if checkpoint_exists "native_${sample}" depth; then
        info "  ${sample}: depth done, skipping"
    else
        local bases
        bases=$(samtools stats -@ 4 -F 0xD04 "${markdup}" \
                | awk -F'\t' '/^SN\tbases mapped \(cigar\)/ {v=$3} END{print v}')
        [[ -n "${bases}" ]] || { echo "ERROR: ${sample}: could not read mapped bases" >&2; return 1; }
        awk -v b="${bases}" -v g="${GENOME_SIZE}" 'BEGIN{printf "%.4f\n", b/g}' > "${MET_DIR}/${sample}.depth.txt"
        log "  ${sample}: $(cat "${MET_DIR}/${sample}.depth.txt")x usable depth"
        create_checkpoint "native_${sample}" depth
    fi
}

phase_align() {
    log "PHASE 1 — align ${#SAMPLES[@]} samples, keep markdup BAMs"
    printf '%s\n' "${SAMPLES[@]}" \
      | PARALLEL_SHELL=/bin/bash parallel --halt soon,fail=1 -j "${MAX_PARALLEL_SAMPLES}" align_sample {}

    # Native vs v2 depth: same reads, same flags; only bwa thread batching differs.
    local cmp="${MET_DIR}/depth_native_vs_v2.tsv"
    printf 'sample\tnative_depth_X\tv2_depth_X\trel_diff\n' > "${cmp}"
    for s in "${SAMPLES[@]}"; do
        local n v; n=$(cat "${MET_DIR}/${s}.depth.txt"); v=$(awk -v s="${s}" '$1==s{print $2}' "${V2_DEPTH}")
        awk -v s="${s}" -v n="${n}" -v v="${v}" 'BEGIN{printf "%s\t%s\t%s\t%.4f\n", s, n, v, (v>0?(n-v)/v:0)}' >> "${cmp}"
    done
    column -t "${cmp}" | while read -r l; do info "    ${l}"; done
    log "PHASE 1 complete"
}

################################################################################
# INTERVALS  (same scheme as generate_all_vcf.sh)
################################################################################

# Writes <name>\t<-L value>\t<length> sorted longest first, to $1
make_intervals() {
    local out=$1
    if [[ -n "${TEST_INTERVAL}" ]]; then
        printf '%s\t%s\t1\n' "$(tr ':-' '__' <<< "${TEST_INTERVAL}")" "${TEST_INTERVAL}" > "${out}"
        return 0
    fi
    awk '$1 !~ /^Chr/ {print $1}' "${REF_GENOME}.fai" > "${NAT_DIR}/scaffolds.list"
    {
        awk '$1 ~ /^Chr/ {print $1"\t"$1"\t"$2}' "${REF_GENOME}.fai"
        awk '$1 !~ /^Chr/ {s+=$2} END{print "scaffolds\t'"${NAT_DIR}"'/scaffolds.list\t"s}' "${REF_GENOME}.fai"
    } | sort -t$'\t' -k3,3nr > "${out}"
}

################################################################################
# PHASE 2 — HAPLOTYPECALLER, scattered sample x interval
################################################################################

hc_job() {
    set -euo pipefail
    local sample=$1 name=$2 iv=$3
    local bam="${BAM_DIR}/${sample}.markdup.bam"
    local gvcf="${NAT_DIR}/gvcf/${sample}/${name}.g.vcf.gz"
    if checkpoint_exists "${CK}_${sample}" "hc_${name}"; then return 0; fi
    [[ -f "${bam}.bai" ]] || { echo "ERROR: no BAM for ${sample}" >&2; return 1; }
    check_disk
    mkdir -p "$(dirname "${gvcf}")" "${TMP_DIR}/${CK}_hc_${sample}_${name}"
    gatk --java-options "-Xmx${MEM_HC}" HaplotypeCaller \
        -R "${REF_GENOME}" -I "${bam}" -L "${iv}" -O "${gvcf}" \
        -ERC GVCF \
        --native-pair-hmm-threads "${HC_THREADS}" \
        --min-base-quality-score "${MIN_BQ}" \
        --tmp-dir "${TMP_DIR}/${CK}_hc_${sample}_${name}" \
        &> "${LOG_DIR}/${CK}_${sample}_${name}_hc.log"
    rm -rf "${TMP_DIR}/${CK}_hc_${sample}_${name}"
    create_checkpoint "${CK}_${sample}" "hc_${name}"
    log "  hc ${sample}/${name} done"
}

phase_call() {
    local iv="${NAT_DIR}/intervals.tsv"; make_intervals "${iv}"
    local jobs="${NAT_DIR}/hc_jobs.tsv"
    # Longest intervals first across all samples, so the tail of the run is short.
    while IFS=$'\t' read -r name L len; do
        while read -r s; do printf '%s\t%s\t%s\n' "${s}" "${name}" "${L}"; done < "${CALL_SAMPLES}"
    done < "${iv}" > "${jobs}"
    log "PHASE 2 — HaplotypeCaller: $(wc -l < "${jobs}") jobs, ${HC_JOBS} at a time x ${HC_THREADS} threads"
    PARALLEL_SHELL=/bin/bash parallel --colsep '\t' --halt soon,fail=1 -j "${HC_JOBS}" \
        hc_job {1} {2} {3} :::: "${jobs}"
    log "PHASE 2 complete"
}

################################################################################
# PHASE 3 — JOINT GENOTYPING per interval, then concat
################################################################################

joint_job() {
    set -euo pipefail
    local name=$1 iv=$2
    local combined="${NAT_DIR}/pieces/${name}.g.vcf.gz"
    local out="${NAT_DIR}/pieces/${name}.vcf.gz"
    if checkpoint_exists "${CK}" "joint_${name}"; then return 0; fi
    check_disk
    mkdir -p "${TMP_DIR}/${CK}_joint_${name}"
    local args=()
    while read -r s; do args+=(-V "${NAT_DIR}/gvcf/${s}/${name}.g.vcf.gz"); done < "${CALL_SAMPLES}"
    gatk --java-options "-Xmx${MEM_JOINT}" CombineGVCFs -R "${REF_GENOME}" "${args[@]}" \
        -L "${iv}" -O "${combined}" --tmp-dir "${TMP_DIR}/${CK}_joint_${name}" \
        &> "${LOG_DIR}/${CK}_${name}_combine.log"
    gatk --java-options "-Xmx${MEM_JOINT}" GenotypeGVCFs -R "${REF_GENOME}" -V "${combined}" \
        -L "${iv}" -O "${out}" --tmp-dir "${TMP_DIR}/${CK}_joint_${name}" \
        &> "${LOG_DIR}/${CK}_${name}_genotype.log"
    rm -f "${combined}" "${combined}.tbi"
    rm -rf "${TMP_DIR}/${CK}_joint_${name}"
    create_checkpoint "${CK}" "joint_${name}"
    log "  joint ${name} done"
}

phase_joint() {
    local iv="${NAT_DIR}/intervals.tsv"; [[ -f "${iv}" ]] || make_intervals "${iv}"
    local raw="${NAT_DIR}/raw/native13.raw.vcf.gz"
    log "PHASE 3 — joint genotyping, $(wc -l < "${iv}") intervals"
    cut -f1,2 "${iv}" | PARALLEL_SHELL=/bin/bash parallel --colsep '\t' --halt soon,fail=1 \
        -j "${JOINT_JOBS}" joint_job {1} {2}

    # Concatenate in reference (.fai) order, not the longest-first job order.
    local pieces=()
    if [[ -n "${TEST_INTERVAL}" ]]; then
        pieces=("${NAT_DIR}/pieces/$(cut -f1 "${iv}").vcf.gz")
    else
        while read -r c; do pieces+=("${NAT_DIR}/pieces/${c}.vcf.gz"); done \
            < <(awk '$1 ~ /^Chr/ {print $1}' "${REF_GENOME}.fai")
        pieces+=("${NAT_DIR}/pieces/scaffolds.vcf.gz")
    fi
    bcftools concat --threads 8 -O z -o "${raw}" "${pieces[@]}"
    tabix -f -p vcf "${raw}"
    log "  $(bcftools view -H "${raw}" | wc -l) variant records, $(bcftools query -l "${raw}" | wc -l) samples -> ${raw}"
    log "PHASE 3 complete"
}

export -f align_sample hc_job joint_job _emit log info warn error checkpoint_exists \
          create_checkpoint disk_pct check_disk
export READS_DIR REF_GENOME WORK_DIR GENOME_SIZE LOG_DIR TMP_DIR CHECKPOINT_DIR MAIN_LOG \
       NAT_DIR CK BAM_DIR MET_DIR CALL_SAMPLES THREADS_PER_SAMPLE MEM_PER_SAMPLE HC_THREADS \
       MEM_HC MEM_JOINT MIN_BQ DISK_WARN_PCT DISK_ABORT_PCT FORCE_RERUN RED GREEN YELLOW BLUE NC

################################################################################
# MAIN
################################################################################

main() {
    preflight
    case "${PHASE}" in
        all)   phase_align; phase_call; phase_joint ;;
        align) phase_align ;;
        call)  phase_call ;;
        joint) phase_joint ;;
        *)     error "unknown phase '${PHASE}' (use: all|align|call|joint)" ;;
    esac
    log "native_depth.sh finished (phase ${PHASE})"
}

main
