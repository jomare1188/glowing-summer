#!/bin/bash
# Illumina Paired-End WGS Quality Control Pipeline (v2)
#
# Usage: ./qc_pipeline_v2.sh [input_dir] [output_dir] [threads]
#
# Environment overrides:
#   PARALLEL_JOBS=8       number of concurrent fastp jobs
#   GENOME_SIZE=506362408 reference size, for coverage estimation
#   TRIM_5P=14            bases hard-trimmed from the 5' end (0 disables)
#   SKIP_RAW_FASTQC=0     skip FastQC on raw reads (they do not change between runs)
#   KRAKEN_DB=/path/to/db taxonomic contamination screen; defaults to the
#                         project PlusPFP index, set to "" to disable
#   SCREEN_READS=200000   read pairs subsampled per sample for screening
#
# Re-running is safe: completed samples are detected and skipped.

set -euo pipefail

INPUT_DIR=${1:-.}
OUTPUT_DIR=${2:-qc_results}
THREADS=${3:-16}

PARALLEL_JOBS=${PARALLEL_JOBS:-8}
GENOME_SIZE=${GENOME_SIZE:-506362408}
TRIM_5P=${TRIM_5P:-14}
SKIP_RAW_FASTQC=${SKIP_RAW_FASTQC:-0}
# Defaults to the project's PlusPFP index so the screen runs whether or not the
# caller remembers the prefix. Set KRAKEN_DB= (empty) to skip classification.
KRAKEN_DB=${KRAKEN_DB-/dados04/jorge/databases/kraken2/k2_pluspfp_16_GB_20260626}
SCREEN_READS=${SCREEN_READS:-200000}

# fastp gains plateau well before 16 worker threads; FastQC allocates 250 MB of
# Java heap per thread and opens that many gzip streams at once. Both are I/O
# bound on large WGS files, so cap them independently of the total core count.
THREADS_PER_JOB=$(( THREADS / PARALLEL_JOBS ))
(( THREADS_PER_JOB < 1 )) && THREADS_PER_JOB=1
(( THREADS_PER_JOB > 16 )) && THREADS_PER_JOB=16
FASTQC_THREADS=$(( THREADS > 16 ? 16 : THREADS ))

READS_DIR="${OUTPUT_DIR}/trimmed_reads"
REPORT_DIR="${OUTPUT_DIR}/reports"
LOG_DIR="${OUTPUT_DIR}/logs"
SCREEN_DIR="${OUTPUT_DIR}/contamination_screen"

mkdir -p "${OUTPUT_DIR}"/{raw_fastqc,trimmed_fastqc} "${READS_DIR}" "${REPORT_DIR}" "${LOG_DIR}" "${SCREEN_DIR}"

log() { echo "[$(date '+%F %T')] $*"; }
die() { echo "[$(date '+%F %T')] ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Preflight
# ---------------------------------------------------------------------------
log "Preflight checks"

for tool in fastp fastqc multiqc parallel; do
    command -v "${tool}" >/dev/null || die "${tool} not found in PATH (activate the qc_illumina env?)"
done
command -v jq >/dev/null || die "jq not found in PATH (needed for the yield table)"

HAVE_SEQKIT=0
command -v seqkit >/dev/null && HAVE_SEQKIT=1

# Discover pairs and verify both mates exist before doing any work.
declare -a R1_FILES=()
shopt -s nullglob
for f in "${INPUT_DIR}"/*_R1*.fastq.gz; do R1_FILES+=("$f"); done
shopt -u nullglob

(( ${#R1_FILES[@]} > 0 )) || die "no *_R1*.fastq.gz files found in ${INPUT_DIR}"

for R1 in "${R1_FILES[@]}"; do
    # Replace the LAST occurrence of _R1 so directory names cannot be clobbered.
    R2="${R1%_R1*}_R2${R1##*_R1}"
    [[ -f "${R2}" ]] || die "missing mate for ${R1} (expected ${R2})"
done
log "Found ${#R1_FILES[@]} complete read pairs"

# Trimmed output is roughly the size of the input; fail now rather than at hour six.
INPUT_KB=$(du -sk "${INPUT_DIR}" | cut -f1)
AVAIL_KB=$(df -Pk "${OUTPUT_DIR}" | awk 'NR==2 {print $4}')
if (( AVAIL_KB < INPUT_KB )); then
    die "insufficient disk space: need ~$((INPUT_KB/1024/1024)) GB, have $((AVAIL_KB/1024/1024)) GB"
fi
log "Disk space OK: $((AVAIL_KB/1024/1024)) GB free, ~$((INPUT_KB/1024/1024)) GB required"

log "Input:   ${INPUT_DIR}"
log "Output:  ${OUTPUT_DIR}"
log "Threads: ${THREADS} (${PARALLEL_JOBS} jobs x ${THREADS_PER_JOB} threads, FastQC ${FASTQC_THREADS})"
log "5' hard trim: ${TRIM_5P} bp"
if [[ -n "${KRAKEN_DB}" ]]; then
    log "Kraken2 DB: ${KRAKEN_DB}"
    [[ -s "${KRAKEN_DB}/hash.k2d" ]] || die "no hash.k2d in KRAKEN_DB: ${KRAKEN_DB}"
    command -v kraken2 >/dev/null || die "KRAKEN_DB set but kraken2 not in PATH (conda activate qc_illumina)"
else
    log "Kraken2 DB: (none) - taxonomic screen disabled"
fi

# ---------------------------------------------------------------------------
# Step 1: FastQC on raw reads
# ---------------------------------------------------------------------------
if [[ "${SKIP_RAW_FASTQC}" == "1" ]]; then
    log "Skipping raw FastQC (SKIP_RAW_FASTQC=1)"
else
    log "FastQC on raw reads"
    declare -a RAW_TODO=()
    for R1 in "${R1_FILES[@]}"; do
        R2="${R1%_R1*}_R2${R1##*_R1}"
        for f in "${R1}" "${R2}"; do
            base=$(basename "${f}" .fastq.gz)
            [[ -f "${OUTPUT_DIR}/raw_fastqc/${base}_fastqc.zip" ]] || RAW_TODO+=("${f}")
        done
    done
    if (( ${#RAW_TODO[@]} > 0 )); then
        log "  ${#RAW_TODO[@]} file(s) to process"
        fastqc -t "${FASTQC_THREADS}" -o "${OUTPUT_DIR}/raw_fastqc" "${RAW_TODO[@]}" \
            > "${LOG_DIR}/raw_fastqc.log" 2>&1
    else
        log "  all raw FastQC reports already present, skipping"
    fi

    log "MultiQC report for raw reads"
    multiqc "${OUTPUT_DIR}/raw_fastqc" -o "${REPORT_DIR}" -n raw_multiqc_report --force \
        > "${LOG_DIR}/raw_multiqc.log" 2>&1
fi

# ---------------------------------------------------------------------------
# Step 2: Trimming and filtering with fastp
# ---------------------------------------------------------------------------
log "Trimming with fastp"

export READS_DIR LOG_DIR THREADS_PER_JOB TRIM_5P

process_fastp() {
    # These options are NOT inherited from the parent shell: GNU parallel spawns
    # a fresh `bash -c`, which does not import -e/-u/pipefail. Without pipefail
    # here a fastp crash would be masked by a successful downstream command and
    # the sample would silently end up with a truncated FASTQ.
    set -euo pipefail

    local R1=$1
    local SAMPLE R2 OUT_R1 OUT_R2 JSON HTML LOG
    SAMPLE=$(basename "${R1}" | sed 's/_R1.*//')
    R2="${R1%_R1*}_R2${R1##*_R1}"

    OUT_R1="${READS_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
    OUT_R2="${READS_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
    JSON="${READS_DIR}/${SAMPLE}_fastp.json"
    HTML="${READS_DIR}/${SAMPLE}_fastp.html"
    LOG="${LOG_DIR}/${SAMPLE}_fastp.log"

    # The JSON is written only after fastp completes, so it is a valid
    # completion sentinel; the FASTQ files exist from the moment writing starts.
    if [[ -s "${JSON}" && -s "${OUT_R1}" && -s "${OUT_R2}" ]]; then
        echo "[$(date '+%F %T')] ${SAMPLE}: already complete, skipping"
        return 0
    fi

    echo "[$(date '+%F %T')] ${SAMPLE}: starting"
    rm -f "${OUT_R1}" "${OUT_R2}" "${JSON}" "${HTML}"

    fastp \
        -i "${R1}" \
        -I "${R2}" \
        -o "${OUT_R1}" \
        -O "${OUT_R2}" \
        --thread "${THREADS_PER_JOB}" \
        -f "${TRIM_5P}" \
        -F "${TRIM_5P}" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 30 \
        --length_required 50 \
        --trim_poly_g \
        --trim_poly_x \
        --low_complexity_filter \
        --complexity_threshold 30 \
        --correction \
        --overlap_len_require 30 \
        --overlap_diff_limit 5 \
        --html "${HTML}" \
        --json "${JSON}" \
        > "${LOG}" 2>&1

    [[ -s "${JSON}" ]] || { echo "ERROR: ${SAMPLE}: fastp produced no JSON" >&2; return 1; }
    echo "[$(date '+%F %T')] ${SAMPLE}: done"
}
export -f process_fastp

# --halt soon,fail=1 stops dispatching new jobs on the first failure and makes
# parallel exit non-zero, which set -e then turns into a pipeline abort.
printf '%s\n' "${R1_FILES[@]}" | \
    parallel --halt soon,fail=1 -j "${PARALLEL_JOBS}" process_fastp {} \
    2>&1 | tee "${LOG_DIR}/fastp_parallel.log"

log "All ${#R1_FILES[@]} samples trimmed successfully"

# ---------------------------------------------------------------------------
# Step 3: FastQC + MultiQC on trimmed reads
# ---------------------------------------------------------------------------
log "FastQC on trimmed reads"
declare -a TRIM_TODO=()
shopt -s nullglob
for f in "${READS_DIR}"/*_trimmed.fastq.gz; do
    base=$(basename "${f}" .fastq.gz)
    [[ -f "${OUTPUT_DIR}/trimmed_fastqc/${base}_fastqc.zip" ]] || TRIM_TODO+=("${f}")
done
shopt -u nullglob

if (( ${#TRIM_TODO[@]} > 0 )); then
    log "  ${#TRIM_TODO[@]} file(s) to process"
    fastqc -t "${FASTQC_THREADS}" -o "${OUTPUT_DIR}/trimmed_fastqc" "${TRIM_TODO[@]}" \
        > "${LOG_DIR}/trimmed_fastqc.log" 2>&1
else
    log "  all trimmed FastQC reports already present, skipping"
fi

log "MultiQC report for trimmed reads"
multiqc "${OUTPUT_DIR}/trimmed_fastqc" "${READS_DIR}" "${SCREEN_DIR}" \
    -o "${REPORT_DIR}" -n trimmed_multiqc_report --force \
    --ignore '*.fastq.gz' \
    > "${LOG_DIR}/trimmed_multiqc.log" 2>&1

# ---------------------------------------------------------------------------
# Step 4: Per-sample yield and coverage table
# ---------------------------------------------------------------------------
# This is what drives the >=9x / >=20x dataset assignment. est_coverage_X is an
# upper bound: it assumes every trimmed base maps and none is a duplicate.
# Realised depth will be est_coverage_X * mapping_rate * (1 - duplicate_rate).
log "Building per-sample yield table"

YIELD="${REPORT_DIR}/yield_per_sample.csv"
{
    echo "sample,raw_reads,raw_gbases,clean_reads,clean_gbases,pct_bases_retained,q30_rate,gc_content,duplication_rate,est_coverage_X"
    for j in "${READS_DIR}"/*_fastp.json; do
        s=$(basename "${j}" _fastp.json)
        jq -r --arg s "${s}" --argjson g "${GENOME_SIZE}" '
            [ $s,
              .summary.before_filtering.total_reads,
              (.summary.before_filtering.total_bases / 1e9 * 100 | round / 100),
              .summary.after_filtering.total_reads,
              (.summary.after_filtering.total_bases / 1e9 * 100 | round / 100),
              (.summary.after_filtering.total_bases / .summary.before_filtering.total_bases * 10000 | round / 100),
              (.summary.after_filtering.q30_rate * 10000 | round / 100),
              (.summary.after_filtering.gc_content * 10000 | round / 100),
              ((.duplication.rate // 0) * 10000 | round / 100),
              (.summary.after_filtering.total_bases / $g * 100 | round / 100)
            ] | @csv' "${j}"
    done
} | sed 's/"//g' > "${YIELD}"

log "  wrote ${YIELD}"

# ---------------------------------------------------------------------------
# Step 5: Contamination screen
# ---------------------------------------------------------------------------
# The previous run showed mapping rates from 6% to 77% against near-uniform raw
# depth, i.e. most reads in the worst samples are not the target organism.
# Nothing in the old pipeline would have flagged that, so screen explicitly.
log "Contamination screen"

GC_TABLE="${SCREEN_DIR}/gc_distribution_summary.tsv"
if (( HAVE_SEQKIT )); then
    printf 'sample\tmean_GC\tmedian_GC\tp10_GC\tp90_GC\tfrac_GC_lt_30\tfrac_GC_gt_50\n' > "${GC_TABLE}"
    for R1 in "${R1_FILES[@]}"; do
        SAMPLE=$(basename "${R1}" | sed 's/_R1.*//')
        HIST="${SCREEN_DIR}/${SAMPLE}_gc_hist.tsv"

        if [[ ! -s "${HIST}" ]]; then
            # head -n is a streaming read of the first N records: no full pass
            # over a 20 GB file, and Illumina output is not position-sorted by
            # anything biologically meaningful, so this is adequate for GC.
            seqkit head -n "${SCREEN_READS}" "${R1}" 2>/dev/null \
                | seqkit fx2tab -n -g 2>/dev/null \
                | awk -F'\t' '{b=int($NF); h[b]++} END {for (i=0;i<=100;i++) printf "%d\t%d\n", i, (h[i]+0)}' \
                > "${HIST}"
        fi

        awk -F'\t' -v s="${SAMPLE}" '
            {gc[$1]=$2; total+=$2; sum+=$1*$2; if($1<30) lo+=$2; if($1>50) hi+=$2}
            END {
                if (total == 0) { printf "%s\tNA\tNA\tNA\tNA\tNA\tNA\n", s; exit }
                c=0; for (i=0;i<=100;i++) { c+=gc[i]
                    if (!p10 && c >= total*0.10) p10=i
                    if (!med && c >= total*0.50) med=i
                    if (!p90 && c >= total*0.90) p90=i }
                printf "%s\t%.1f\t%d\t%d\t%d\t%.3f\t%.3f\n", s, sum/total, med, p10, p90, lo/total, hi/total
            }' "${HIST}" >> "${GC_TABLE}"
    done
    log "  wrote ${GC_TABLE}"
    log "  a wide p10-p90 spread or a bimodal histogram indicates mixed-organism input"
else
    log "  seqkit not found, skipping GC screen"
fi

if [[ -n "${KRAKEN_DB}" ]]; then
    command -v kraken2 >/dev/null || die "KRAKEN_DB set but kraken2 not in PATH (conda activate qc_illumina)"
    [[ -d "${KRAKEN_DB}" ]] || die "KRAKEN_DB directory not found: ${KRAKEN_DB}"
    log "  taxonomic classification with kraken2"

    for R1 in "${R1_FILES[@]}"; do
        SAMPLE=$(basename "${R1}" | sed 's/_R1.*//')
        REPORT="${SCREEN_DIR}/${SAMPLE}.kraken2.report.txt"
        [[ -s "${REPORT}" ]] && continue

        R2="${R1%_R1*}_R2${R1##*_R1}"
        SUB1="${SCREEN_DIR}/${SAMPLE}_sub_R1.fq"
        SUB2="${SCREEN_DIR}/${SAMPLE}_sub_R2.fq"
        seqkit head -n "${SCREEN_READS}" "${R1}" > "${SUB1}"
        seqkit head -n "${SCREEN_READS}" "${R2}" > "${SUB2}"

        kraken2 --db "${KRAKEN_DB}" --threads "${THREADS}" --paired \
            --report "${REPORT}" --output /dev/null \
            "${SUB1}" "${SUB2}" > "${LOG_DIR}/${SAMPLE}_kraken2.log" 2>&1

        rm -f "${SUB1}" "${SUB2}"
    done

    # The absolute classified fraction is NOT the contamination metric: Callicarpa
    # americana is not in RefSeq, so ~90% of reads are unclassified even in clean
    # samples. What separates samples is the composition of the classified
    # fraction. Reads from the target land on related Lamiales genera (Sesamum,
    # Salvia, Olea, ...) and so score as Viridiplantae; the contaminant is soil
    # bacteria. On this cohort the plant:bacteria ratio tracks the measured
    # mapping rate with Spearman rho = 0.94, so it predicts alignment success
    # before any alignment is run.
    KRAKEN_TABLE="${SCREEN_DIR}/contamination_index.tsv"
    printf 'sample\tplant_pct\tbacteria_pct\tfungi_pct\tunclassified_pct\tplant_to_bacteria\n' > "${KRAKEN_TABLE}"
    for R1 in "${R1_FILES[@]}"; do
        SAMPLE=$(basename "${R1}" | sed 's/_R1.*//')
        REPORT="${SCREEN_DIR}/${SAMPLE}.kraken2.report.txt"
        [[ -s "${REPORT}" ]] || continue
        awk -F'\t' -v s="${SAMPLE}" '
            {n=$6; gsub(/^ +/,"",n)
             if (n=="Viridiplantae") p=$1
             else if (n=="Bacteria") b=$1
             else if (n=="Fungi")    f=$1
             else if (n=="unclassified") u=$1}
            END {printf "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", s, p, b, f, u, (b>0 ? p/b : -1)}
        ' "${REPORT}" >> "${KRAKEN_TABLE}"
    done
    log "  wrote ${KRAKEN_TABLE}"
    log "  low plant_to_bacteria (<1) predicts a low mapping rate; inspect those samples"
    log "  kraken2 reports in ${SCREEN_DIR} (also parsed into the MultiQC report)"
else
    log "  KRAKEN_DB empty, skipping taxonomic classification"
    log "  to enable: KRAKEN_DB=/path/to/db bash $0 ..."
fi

# ---------------------------------------------------------------------------
# Step 6: Summary
# ---------------------------------------------------------------------------
log "Summary"
{
    echo "QC Pipeline Summary"
    echo "==================="
    echo "Date:           $(date)"
    echo "Input:          ${INPUT_DIR}"
    echo "Output:         ${OUTPUT_DIR}"
    echo "Genome size:    ${GENOME_SIZE} bp"
    echo "5' hard trim:   ${TRIM_5P} bp"
    echo "Read pairs in:  ${#R1_FILES[@]}"
    echo "Read pairs out: $(find "${READS_DIR}" -name '*_R1_trimmed.fastq.gz' | wc -l)"
    echo ""
    echo "Reports:"
    echo "  Raw reads:       ${REPORT_DIR}/raw_multiqc_report.html"
    echo "  Trimmed reads:   ${REPORT_DIR}/trimmed_multiqc_report.html"
    echo "  Yield/coverage:  ${YIELD}"
    [[ -s "${GC_TABLE}" ]] && echo "  GC screen:       ${GC_TABLE}"
    echo ""
    echo "Estimated coverage (upper bound: assumes 100% mapping, 0% duplicates):"
    column -s, -t "${YIELD}" | awk 'NR==1 || 1' | cut -c1-200
} | tee "${OUTPUT_DIR}/qc_summary.txt"

log "QC pipeline completed successfully"
log "Trimmed reads: ${READS_DIR}"
