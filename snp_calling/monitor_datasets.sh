#!/bin/bash
################################################################################
# MONITOR for snp_datasets.sh
#
# Read-only: inspects checkpoints, logs and output files. It never writes into
# the results tree, so it is safe to run while the pipeline is going.
#
# Usage:  bash monitor_datasets.sh [-w] [-n SECONDS]
#           -w            watch mode, refresh until interrupted
#           -n SECONDS    refresh interval in watch mode (default 30)
################################################################################

set -uo pipefail

PROJECT_DIR="/dados04/jorge/CALLICARPA_2"
WORK_DIR="${PROJECT_DIR}/snp_calling"
OUT_DIR="${WORK_DIR}/results"
LOG_DIR="${WORK_DIR}/logs"
CHECKPOINT_DIR="${WORK_DIR}/checkpoints"
MAIN_LOG="${LOG_DIR}/pipeline_main.log"
READS_DIR="${PROJECT_DIR}/qc/qc_results/trimmed_reads"

DATASETS=("d09:9" "d20:20")

WATCH=false
INTERVAL=30
while getopts "wn:" opt; do
    case $opt in
        w) WATCH=true ;;
        n) INTERVAL=$OPTARG ;;
        *) echo "usage: $0 [-w] [-n SECONDS]" >&2; exit 1 ;;
    esac
done

BOLD='\033[1m'; DIM='\033[2m'; RED='\033[0;31m'; GREEN='\033[0;32m'
YELLOW='\033[1;33m'; BLUE='\033[0;34m'; CYAN='\033[0;36m'; NC='\033[0m'

hr() { printf "${DIM}%s${NC}\n" "$(printf '─%.0s' {1..78})"; }

human() {
    local kb=$1
    awk -v k="$kb" 'BEGIN{
        split("KB MB GB TB", u, " "); i=1
        while (k >= 1024 && i < 4) { k/=1024; i++ }
        printf "%.1f %s", k, u[i]
    }'
}

bar() {
    local done=$1 total=$2 width=${3:-28}
    (( total == 0 )) && { printf "[%${width}s] n/a" ""; return; }
    local filled=$(( done * width / total ))
    local pct=$(( done * 100 / total ))
    printf "["
    printf "%${filled}s" "" | tr ' ' '#'
    printf "%$(( width - filled ))s" "" | tr ' ' '.'
    printf "] %3d%% (%d/%d)" "$pct" "$done" "$total"
}

count_cp() { ls "${CHECKPOINT_DIR}"/*."$1".done 2>/dev/null | wc -l | tr -d ' '; }

n_samples() {
    ls "${READS_DIR}"/*_R1_trimmed.fastq.gz 2>/dev/null | wc -l | tr -d ' '
}

################################################################################

show_header() {
    local state="${RED}not running${NC}"
    local pids
    pids=$(pgrep -f "snp_datasets.sh" | tr '\n' ' ')
    [[ -n "${pids// /}" ]] && state="${GREEN}running${NC} ${DIM}(pid ${pids% })${NC}"

    hr
    printf "${BOLD} SNP DATASET PIPELINE${NC}   %s   ${DIM}%s${NC}\n" "$(echo -e "$state")" "$(date '+%F %T')"
    hr
}

show_storage() {
    printf "${BOLD}STORAGE${NC}\n"
    local pct avail
    pct=$(df -P "${WORK_DIR}" 2>/dev/null | awk 'NR==2{gsub(/%/,"",$5);print $5}')
    avail=$(df -Ph "${WORK_DIR}" 2>/dev/null | awk 'NR==2{print $4}')
    local colour=$GREEN
    (( pct >= 85 )) && colour=$YELLOW
    (( pct >= 95 )) && colour=$RED
    printf "  filesystem   ${colour}%s%%${NC} used, %s free\n" "$pct" "$avail"

    for d in alignment markdup downsampled metrics; do
        if [[ -d "${OUT_DIR}/${d}" ]]; then
            local kb; kb=$(du -sk "${OUT_DIR}/${d}" 2>/dev/null | cut -f1)
            printf "  %-12s %10s\n" "$d" "$(human "${kb:-0}")"
        fi
    done
    for t in "${DATASETS[@]}"; do
        local tag="${t%%:*}"
        [[ -d "${OUT_DIR}/${tag}" ]] || continue
        local kb; kb=$(du -sk "${OUT_DIR}/${tag}" 2>/dev/null | cut -f1)
        printf "  %-12s %10s\n" "$tag" "$(human "${kb:-0}")"
    done
    echo
}

show_phase1() {
    local total; total=$(n_samples)
    printf "${BOLD}PHASE 1 — alignment / markdup / depth${NC}   ${DIM}%s samples${NC}\n" "$total"
    for step in alignment markdup depth; do
        printf "  %-10s %b\n" "$step" "$(bar "$(count_cp "$step")" "$total")"
    done

    # Anything currently mid-flight has a checkpoint for the previous step only.
    local running
    running=$(pgrep -af "bwa mem|MarkDuplicates|samtools stats" 2>/dev/null \
              | grep -oE "ID:[0-9]+|markdup/[0-9]+|metrics/[0-9]+" | grep -oE "[0-9]+$" \
              | sort -u | tr '\n' ' ')
    [[ -n "${running// /}" ]] && printf "  ${CYAN}active${NC}     %s\n" "$running"
    echo
}

show_depths() {
    local table="${OUT_DIR}/metrics/depth_summary.tsv"
    [[ -f "$table" ]] || return 0
    printf "${BOLD}MEASURED DEPTH AND DATASET MEMBERSHIP${NC}\n"
    column -t "$table" | sed 's/^/  /' \
      | awk -v g="$GREEN" -v r="$DIM" -v n="$NC" \
            'NR==1{print; next} {gsub(/\bIN\b/, g "IN" n); gsub(/\bout\b/, r "out" n); print}'
    echo
}

show_datasets() {
    for t in "${DATASETS[@]}"; do
        local tag target
        tag="${t%%:*}"; target="${t##*:}"
        local list="${OUT_DIR}/metrics/samples_${tag}.txt"
        local n=0
        [[ -f "$list" ]] && n=$(wc -l < "$list" | tr -d ' ')

        printf "${BOLD}DATASET %s${NC} ${DIM}(target %sx, %s samples)${NC}\n" "$tag" "$target" "$n"
        if (( n == 0 )); then
            printf "  ${DIM}not yet assigned${NC}\n\n"
            continue
        fi

        printf "  %-14s %b\n" "downsample" "$(bar "$(count_cp "${tag}_downsample")" "$n")"
        printf "  %-14s %b\n" "gvcf"       "$(bar "$(count_cp "${tag}_gvcf")" "$n")"

        for step in joint_genotyping filtering ld_pruning subset_50K subset_5K roh_input; do
            local mark="${DIM}pending${NC}"
            [[ -f "${CHECKPOINT_DIR}/${tag}.${step}.done" ]] && mark="${GREEN}done${NC}"
            printf "  %-14s %b\n" "$step" "$mark"
        done

        # Output inventory with SNP counts, when the files exist.
        local dd="${OUT_DIR}/${tag}"
        for f in "${dd}/filtered/cohort.snps.final.vcf.gz" \
                 "${dd}/downstream/ld_pruned.snps.vcf.gz" \
                 "${dd}/downstream/ld_pruned_50K.snps.vcf.gz" \
                 "${dd}/downstream/ld_pruned_5K.snps.vcf.gz" \
                 "${dd}/roh/for_roh_dense_maf005.vcf.gz" \
                 "${dd}/roh/for_roh_dense_maf010.vcf.gz"; do
            [[ -f "$f" ]] || continue
            local cnt="?"
            [[ -f "${f}.tbi" || -f "${f}.csi" ]] && cnt=$(bcftools index -n "$f" 2>/dev/null || echo "?")
            printf "    ${DIM}%-42s${NC} %12s SNPs\n" "$(basename "$f")" "$cnt"
        done
        echo
    done
}

show_errors() {
    local errs
    errs=$(grep -h "\[ERROR\]" "${MAIN_LOG}" 2>/dev/null | tail -3)
    if [[ -n "$errs" ]]; then
        printf "${BOLD}${RED}ERRORS${NC}\n"
        echo "$errs" | sed 's/^/  /'
        echo
    fi

    # A GATK or bwa log that ends in an exception is the usual failure mode and
    # does not always reach the main log.
    local bad=()
    shopt -s nullglob
    for l in "${LOG_DIR}"/*.log; do
        [[ "$l" == "${MAIN_LOG}" ]] && continue
        if tail -25 "$l" 2>/dev/null | grep -qE "Exception|ERROR|Traceback|Killed|No space left"; then
            bad+=("$(basename "$l")")
        fi
    done
    shopt -u nullglob
    if (( ${#bad[@]} > 0 )); then
        printf "${BOLD}${YELLOW}LOGS WITH ERROR LINES${NC} ${DIM}(%d)${NC}\n" "${#bad[@]}"
        printf '  %s\n' "${bad[@]:0:6}"
        (( ${#bad[@]} > 6 )) && printf "  ${DIM}... and %d more${NC}\n" "$(( ${#bad[@]} - 6 ))"
        echo
    fi
}

show_recent() {
    [[ -f "${MAIN_LOG}" ]] || return 0
    printf "${BOLD}RECENT ACTIVITY${NC}\n"
    tail -8 "${MAIN_LOG}" | sed 's/^/  /'
    echo
}

render() {
    show_header
    show_storage
    show_phase1
    show_depths
    show_datasets
    show_errors
    show_recent
    hr
}

if $WATCH; then
    trap 'printf "\n"; exit 0' INT
    while true; do
        clear
        render
        printf "${DIM} refreshing every %ss — Ctrl-C to exit${NC}\n" "$INTERVAL"
        sleep "$INTERVAL"
    done
else
    render
fi
