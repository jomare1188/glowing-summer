#!/bin/bash

################################################################################
# SNP Pipeline Monitor - SLURM Compatible
# Monitors pipeline progress via log files and checkpoint files only
# Can run on different machine from the main pipeline
################################################################################

# Configuration (match your pipeline settings)
WORK_DIR="/home/dmpachon/jorge/TATIANA/snp_calling/"
OUT_DIR="${WORK_DIR}/results"
LOG_DIR="${WORK_DIR}/logs_parallel"
CHECKPOINT_DIR="${WORK_DIR}/checkpoints"
READS_DIR="/home/dmpachon/jorge/TATIANA/qc/qc_results/trimmed_reads"
MAIN_LOG="${LOG_DIR}/pipeline_main.log"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m'
BOLD='\033[1m'

################################################################################
# HELPER FUNCTIONS
################################################################################

human_readable_size() {
    local size=$1
    if [[ $size -gt 1073741824 ]]; then
        echo "$(awk "BEGIN {printf \"%.2f\", $size/1073741824}")GB"
    elif [[ $size -gt 1048576 ]]; then
        echo "$(awk "BEGIN {printf \"%.2f\", $size/1048576}")MB"
    elif [[ $size -gt 1024 ]]; then
        echo "$(awk "BEGIN {printf \"%.2f\", $size/1024}")KB"
    else
        echo "${size}B"
    fi
}

get_dir_size() {
    local dir=$1
    if [[ -d "$dir" ]]; then
        du -sb "$dir" 2>/dev/null | awk '{print $1}'
    else
        echo "0"
    fi
}

count_files() {
    local dir=$1
    local pattern=$2
    if [[ -d "$dir" ]]; then
        find "$dir" -name "$pattern" 2>/dev/null | wc -l
    else
        echo "0"
    fi
}

get_file_age() {
    local file=$1
    if [[ -f "$file" ]]; then
        local now=$(date +%s)
        local file_time=$(stat -c %Y "$file" 2>/dev/null || stat -f %m "$file" 2>/dev/null)
        local age=$((now - file_time))
        
        if [[ $age -lt 60 ]]; then
            echo "${age}s ago"
        elif [[ $age -lt 3600 ]]; then
            echo "$((age / 60))m ago"
        elif [[ $age -lt 86400 ]]; then
            echo "$((age / 3600))h ago"
        else
            echo "$((age / 86400))d ago"
        fi
    else
        echo "N/A"
    fi
}

################################################################################
# LOG PARSING FUNCTIONS
################################################################################

get_pipeline_status() {
    if [[ ! -f "$MAIN_LOG" ]]; then
        echo "NOT_STARTED"
        return
    fi
    
    local log_age_seconds=$(( $(date +%s) - $(stat -c %Y "$MAIN_LOG" 2>/dev/null || stat -f %m "$MAIN_LOG" 2>/dev/null) ))
    
    # Check for completion
    if grep -q "Pipeline complete" "$MAIN_LOG" 2>/dev/null; then
        echo "COMPLETED"
    elif grep -q "ERROR" "$MAIN_LOG" 2>/dev/null && [[ $log_age_seconds -gt 300 ]]; then
        echo "FAILED"
    elif [[ $log_age_seconds -lt 300 ]]; then
        echo "RUNNING"
    else
        echo "STALLED"
    fi
}

get_current_step() {
    if [[ ! -f "$MAIN_LOG" ]]; then
        echo "Pipeline not started"
        return
    fi
    
    local last_entry=$(tac "$MAIN_LOG" 2>/dev/null | grep -m1 "^\[" | head -1)
    
    if [[ -n "$last_entry" ]]; then
        echo "$last_entry"
    else
        echo "No recent activity"
    fi
}

parse_slurm_status() {
    if command -v squeue &>/dev/null; then
        local job_info
        job_info=$(squeue -u "$USER" -n snp_calling 2>/dev/null)

        if [[ -n "$job_info" ]] && [[ $(wc -l <<< "$job_info") -gt 1 ]]; then
            local running
            local pending
            local total
            local job_id
            local status_msg="SLURM: "

            running=$(grep -c " R " <<< "$job_info")
            pending=$(grep -c " PD " <<< "$job_info")
            total=$(tail -n +2 <<< "$job_info" | wc -l)

            (( running > 0 )) && status_msg+="${running} running "
            (( pending > 0 )) && status_msg+="${pending} pending "

            job_id=$(tail -n +2 <<< "$job_info" | awk '{print $1}' | head -1)
            [[ -n "$job_id" ]] && status_msg+="(Job ID: $job_id)"

            echo "$status_msg"
        else
            echo "SLURM: No jobs found"
        fi
    else
        echo "SLURM: squeue not available"
    fi
}

################################################################################
# DISPLAY FUNCTIONS
################################################################################

print_header() {
    clear
    echo -e "${BOLD}${CYAN}╔════════════════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BOLD}${CYAN}║          SNP PIPELINE MONITOR - SLURM Compatible                          ║${NC}"
    echo -e "${BOLD}${CYAN}╚════════════════════════════════════════════════════════════════════════════╝${NC}"
    echo -e "${BOLD}Monitoring started: $(date '+%Y-%m-%d %H:%M:%S')${NC}"
    echo ""
}

show_pipeline_status() {
    echo -e "${BOLD}${GREEN}[1] PIPELINE STATUS${NC}"
    echo -e "${BOLD}────────────────────────────────────────────────────────────────────────────${NC}"
    
    local status=$(get_pipeline_status)
    local slurm_status=$(parse_slurm_status)
    
    case $status in
        "RUNNING")
            echo -e "  Pipeline: ${GREEN}● RUNNING${NC}"
            ;;
        "COMPLETED")
            echo -e "  Pipeline: ${GREEN}✓ COMPLETED${NC}"
            ;;
        "FAILED")
            echo -e "  Pipeline: ${RED}✗ FAILED${NC}"
            ;;
        "STALLED")
            echo -e "  Pipeline: ${YELLOW}⚠ STALLED (no recent activity)${NC}"
            ;;
        *)
            echo -e "  Pipeline: ${BLUE}○ NOT STARTED${NC}"
            ;;
    esac
    
    echo -e "  $slurm_status"
    
    if [[ -f "$MAIN_LOG" ]]; then
        echo -e "  Main log: $(basename "$MAIN_LOG") (updated $(get_file_age "$MAIN_LOG"))"
    fi
    
    echo ""
}

show_storage_status() {
    echo -e "${BOLD}${MAGENTA}[2] STORAGE STATUS${NC}"
    echo -e "${BOLD}────────────────────────────────────────────────────────────────────────────${NC}"
    
    if [[ ! -d "$WORK_DIR" ]]; then
        echo -e "  ${RED}✗ Work directory not accessible: ${WORK_DIR}${NC}"
        echo ""
        return
    fi
    
    local disk_info=$(df -h "${WORK_DIR}" 2>/dev/null | tail -n 1)
    if [[ -n "$disk_info" ]]; then
        local disk_use=$(echo "$disk_info" | awk '{print $5}' | head -n 1 | sed 's/%//' | tr -d '[:space:]')
        local disk_avail=$(echo "$disk_info" | awk '{print $4}' | head -n 1 | tr -d '[:space:]')
        local disk_total=$(echo "$disk_info" | awk '{print $2}' | head -n 1 | tr -d '[:space:]')
        
        # Validate disk_use is a number
        if [[ "$disk_use" =~ ^[0-9]+$ ]]; then
            if [[ $disk_use -gt 90 ]]; then
                echo -e "  ${RED}⚠ Disk Usage:${NC} ${RED}${disk_use}%${NC} (${disk_avail} of ${disk_total} available)"
            elif [[ $disk_use -gt 80 ]]; then
                echo -e "  ${YELLOW}⚠ Disk Usage:${NC} ${YELLOW}${disk_use}%${NC} (${disk_avail} of ${disk_total} available)"
            else
                echo -e "  ${GREEN}✓ Disk Usage:${NC} ${disk_use}% (${disk_avail} of ${disk_total} available)"
            fi
        else
            echo -e "  Disk Usage: Unable to parse (${disk_avail} available)"
        fi
    fi
    
    local alignment_size=$(get_dir_size "${OUT_DIR}/alignment")
    local markdup_size=$(get_dir_size "${OUT_DIR}/markdup")
    local variants_size=$(get_dir_size "${OUT_DIR}/variants")
    local logs_size=$(get_dir_size "${LOG_DIR}")
    local total_size=$((alignment_size + markdup_size + variants_size + logs_size))
    
    echo ""
    echo -e "  ${BOLD}Storage breakdown:${NC}"
    echo -e "    Alignment BAMs:     $(human_readable_size $alignment_size)"
    echo -e "    Duplicate-marked:   $(human_readable_size $markdup_size)"
    echo -e "    Variants (VCF):     $(human_readable_size $variants_size)"
    echo -e "    Logs:               $(human_readable_size $logs_size)"
    echo -e "    ${BOLD}Total:              $(human_readable_size $total_size)${NC}"
    
    echo ""
}

show_sample_progress() {
    echo -e "${BOLD}${BLUE}[3] SAMPLE PROGRESS${NC}"
    echo -e "${BOLD}────────────────────────────────────────────────────────────────────────────${NC}"
    
    local total_samples=0
    shopt -s nullglob
    for r1 in "${READS_DIR}"/*_R1*.fastq.gz "${READS_DIR}"/*_1.fastq.gz; do
        if [[ -f "$r1" ]]; then
            ((total_samples++))
        fi
    done
    shopt -u nullglob
    
    if [[ $total_samples -eq 0 ]]; then
        echo -e "  ${YELLOW}No samples found in reads directory${NC}"
        echo ""
        return
    fi
    
    local aligned=0
    local markdup=0
    local called=0
    
    if [[ -d "${CHECKPOINT_DIR}" ]]; then
        aligned=$(ls "${CHECKPOINT_DIR}"/*.alignment.done 2>/dev/null | wc -l | tr -d '[:space:]')
        markdup=$(ls "${CHECKPOINT_DIR}"/*.markdup.done 2>/dev/null | wc -l | tr -d '[:space:]')
        called=$(ls "${CHECKPOINT_DIR}"/*.variant_calling.done 2>/dev/null | wc -l | tr -d '[:space:]')
        
        # Ensure they are numbers
        [[ ! "$aligned" =~ ^[0-9]+$ ]] && aligned=0
        [[ ! "$markdup" =~ ^[0-9]+$ ]] && markdup=0
        [[ ! "$called" =~ ^[0-9]+$ ]] && called=0
    fi
    
    local align_pct=0
    local markdup_pct=0
    local called_pct=0
    
    if [[ $total_samples -gt 0 ]]; then
        align_pct=$(awk "BEGIN {printf \"%.1f\", ($aligned/$total_samples)*100}")
        markdup_pct=$(awk "BEGIN {printf \"%.1f\", ($markdup/$total_samples)*100}")
        called_pct=$(awk "BEGIN {printf \"%.1f\", ($called/$total_samples)*100}")
    fi
    
    echo -e "  Total samples: ${BOLD}${total_samples}${NC}"
    echo ""
    echo -e "  ${CYAN}Step 1 - Alignment:${NC}        ${aligned}/${total_samples} (${align_pct}%)"
    draw_progress_bar $aligned $total_samples
    echo ""
    echo -e "  ${CYAN}Step 2 - Mark Duplicates:${NC} ${markdup}/${total_samples} (${markdup_pct}%)"
    draw_progress_bar $markdup $total_samples
    echo ""
    echo -e "  ${CYAN}Step 3 - Variant Calling:${NC} ${called}/${total_samples} (${called_pct}%)"
    draw_progress_bar $called $total_samples
    echo ""
    
    local bam_count=$(count_files "${OUT_DIR}/alignment" "*.sorted.bam" | tr -d '[:space:]')
    local markdup_count=$(count_files "${OUT_DIR}/markdup" "*.markdup.bam" | tr -d '[:space:]')
    local gvcf_count=$(count_files "${OUT_DIR}/variants/raw" "*.g.vcf.gz" | tr -d '[:space:]')
    
    [[ ! "$bam_count" =~ ^[0-9]+$ ]] && bam_count=0
    [[ ! "$markdup_count" =~ ^[0-9]+$ ]] && markdup_count=0
    [[ ! "$gvcf_count" =~ ^[0-9]+$ ]] && gvcf_count=0
    
    echo -e "  ${BOLD}Output files generated:${NC}"
    echo -e "    Sorted BAMs:      ${bam_count}"
    echo -e "    Markdup BAMs:     ${markdup_count}"
    echo -e "    GVCFs:            ${gvcf_count}"
    echo ""
}

draw_progress_bar() {
    local completed=$1
    local total=$2
    local width=50
    
    if [[ $total -eq 0 ]]; then
        echo -e "  [${RED}$(printf '%*s' $width | tr ' ' '░')${NC}] 0%"
        return
    fi
    
    local filled=$((completed * width / total))
    local empty=$((width - filled))
    
    local bar=""
    if [[ $filled -gt 0 ]]; then
        bar+="${GREEN}$(printf '%*s' $filled | tr ' ' '█')${NC}"
    fi
    if [[ $empty -gt 0 ]]; then
        bar+="$(printf '%*s' $empty | tr ' ' '░')"
    fi
    
    local pct=$(awk "BEGIN {printf \"%.1f\", ($completed/$total)*100}")
    echo -e "  [${bar}] ${pct}%"
}

show_joint_genotyping() {
    echo -e "${BOLD}${YELLOW}[4] JOINT GENOTYPING & FILTERING${NC}"
    echo -e "${BOLD}────────────────────────────────────────────────────────────────────────────${NC}"
    
    local joint_done=false
    local filter_done=false
    
    [[ -f "${CHECKPOINT_DIR}/cohort.joint_genotyping.done" ]] && joint_done=true
    [[ -f "${CHECKPOINT_DIR}/cohort.filtering.done" ]] && filter_done=true
    
    local combined_gvcf="${OUT_DIR}/variants/raw/cohort.g.vcf.gz"
    local raw_vcf="${OUT_DIR}/variants/raw/cohort.raw.vcf.gz"
    local filtered_vcf="${OUT_DIR}/variants/filtered/cohort.snps.filtered.vcf.gz"
    local pass_vcf="${OUT_DIR}/variants/filtered/cohort.snps.pass.vcf.gz"
    
    if [[ -f "$combined_gvcf" ]]; then
        local size=$(stat -c%s "$combined_gvcf" 2>/dev/null || stat -f%z "$combined_gvcf" 2>/dev/null || echo "0")
        local age=$(get_file_age "$combined_gvcf")
        echo -e "  ${GREEN}✓${NC} Combined GVCF:    $(human_readable_size $size) ($age)"
    else
        echo -e "  ${YELLOW}○${NC} Combined GVCF:    Not created yet"
    fi
    
    if [[ -f "$raw_vcf" ]]; then
        local size=$(stat -c%s "$raw_vcf" 2>/dev/null || stat -f%z "$raw_vcf" 2>/dev/null || echo "0")
        local age=$(get_file_age "$raw_vcf")
        echo -e "  ${GREEN}✓${NC} Raw VCF:          $(human_readable_size $size) ($age)"
    else
        echo -e "  ${YELLOW}○${NC} Raw VCF:          Not created yet"
    fi
    
    if [[ -f "$filtered_vcf" ]]; then
        local size=$(stat -c%s "$filtered_vcf" 2>/dev/null || stat -f%z "$filtered_vcf" 2>/dev/null || echo "0")
        local age=$(get_file_age "$filtered_vcf")
        echo -e "  ${GREEN}✓${NC} Filtered VCF:     $(human_readable_size $size) ($age)"
    else
        echo -e "  ${YELLOW}○${NC} Filtered VCF:     Not created yet"
    fi
    
    if [[ -f "$pass_vcf" ]]; then
        local size=$(stat -c%s "$pass_vcf" 2>/dev/null || stat -f%z "$pass_vcf" 2>/dev/null || echo "0")
        local age=$(get_file_age "$pass_vcf")
        
        local variant_count="unknown"
        if command -v bcftools &> /dev/null; then
            variant_count=$(bcftools view -H "$pass_vcf" 2>/dev/null | wc -l | tr -d ' ')
        fi
        
        echo -e "  ${GREEN}✓${NC} PASS VCF:         $(human_readable_size $size) (${variant_count} variants, $age)"
    else
        echo -e "  ${YELLOW}○${NC} PASS VCF:         Not created yet"
    fi
    
    echo ""
    if $joint_done; then
        echo -e "  Joint Genotyping: ${GREEN}✓ COMPLETE${NC}"
    else
        echo -e "  Joint Genotyping: ${YELLOW}○ Pending${NC}"
    fi
    
    if $filter_done; then
        echo -e "  Filtering:        ${GREEN}✓ COMPLETE${NC}"
    else
        echo -e "  Filtering:        ${YELLOW}○ Pending${NC}"
    fi
    
    echo ""
}

show_recent_logs() {
    echo -e "${BOLD}${CYAN}[5] RECENT LOG ACTIVITY${NC}"
    echo -e "${BOLD}────────────────────────────────────────────────────────────────────────────${NC}"
    
    if [[ -f "$MAIN_LOG" ]]; then
        echo -e "  ${BOLD}Latest entries from: $(basename "$MAIN_LOG")${NC}"
        echo ""
        
        tail -8 "$MAIN_LOG" 2>/dev/null | while read -r line; do
            echo "    $line"
        done
    else
        echo -e "  ${YELLOW}No log file found${NC}"
    fi
    
    echo ""
    
    if [[ -f "$MAIN_LOG" ]]; then
        local error_count=$(grep -c "ERROR" "$MAIN_LOG" 2>/dev/null || echo "0")
        local warning_count=$(grep -c "WARNING" "$MAIN_LOG" 2>/dev/null || echo "0")
        
        # Remove any whitespace or newlines
        error_count=$(echo "$error_count" | tr -d '[:space:]')
        warning_count=$(echo "$warning_count" | tr -d '[:space:]')
        
        if [[ "$error_count" =~ ^[0-9]+$ ]] && [[ $error_count -gt 0 ]]; then
            echo -e "  ${RED}⚠ Errors found: ${error_count}${NC}"
            echo -e "  ${BOLD}Last error:${NC}"
            grep "ERROR" "$MAIN_LOG" 2>/dev/null | tail -1 | sed 's/^/    /'
            echo ""
        fi
        
        if [[ "$warning_count" =~ ^[0-9]+$ ]] && [[ $warning_count -gt 0 ]]; then
            echo -e "  ${YELLOW}Warnings: ${warning_count}${NC}"
        fi
    fi
    
    echo ""
}

show_sample_details() {
    echo -e "${BOLD}${MAGENTA}[6] SAMPLE DETAILS${NC}"
    echo -e "${BOLD}────────────────────────────────────────────────────────────────────────────${NC}"
    
    echo -e "  ${BOLD}Sample Status (most recent 15):${NC}"
    echo ""
    
    local recent_samples=()
    if [[ -d "${CHECKPOINT_DIR}" ]]; then
        while IFS= read -r checkpoint; do
            if [[ -n "$checkpoint" ]]; then
                local sample=$(basename "$checkpoint" | sed 's/\..*//')
                recent_samples+=("$sample")
            fi
        done < <(ls -t "${CHECKPOINT_DIR}"/*.done 2>/dev/null | head -15)
    fi
    
    if [[ ${#recent_samples[@]} -eq 0 ]]; then
        local count=0
        shopt -s nullglob
        for r1 in "${READS_DIR}"/*_R1*.fastq.gz "${READS_DIR}"/*_1.fastq.gz; do
            if [[ -f "$r1" ]] && [[ $count -lt 15 ]]; then
                local sample=$(basename "$r1" | sed -E 's/_R?1.*\.fastq\.gz$//')
                recent_samples+=("$sample")
                ((count++))
            fi
        done
        shopt -u nullglob
    fi
    
    local shown=0
    for sample in "${recent_samples[@]}"; do
        if [[ $shown -lt 15 ]]; then
            local align_done=" "
            local markdup_done=" "
            local variant_done=" "
            
            [[ -f "${CHECKPOINT_DIR}/${sample}.alignment.done" ]] && align_done="✓"
            [[ -f "${CHECKPOINT_DIR}/${sample}.markdup.done" ]] && markdup_done="✓"
            [[ -f "${CHECKPOINT_DIR}/${sample}.variant_calling.done" ]] && variant_done="✓"
            
            local latest_time="N/A"
            local latest_checkpoint=$(ls -t "${CHECKPOINT_DIR}/${sample}".*.done 2>/dev/null | head -1)
            if [[ -n "$latest_checkpoint" ]]; then
                latest_time=$(get_file_age "$latest_checkpoint")
            fi
            
            if [[ "$variant_done" == "✓" ]]; then
                echo -e "    ${GREEN}${sample}${NC}: [${align_done}] [${markdup_done}] [${variant_done}] - ${GREEN}COMPLETE${NC} ($latest_time)"
            elif [[ "$markdup_done" == "✓" ]]; then
                echo -e "    ${YELLOW}${sample}${NC}: [${align_done}] [${markdup_done}] [${variant_done}] - ${YELLOW}IN PROGRESS${NC} ($latest_time)"
            elif [[ "$align_done" == "✓" ]]; then
                echo -e "    ${CYAN}${sample}${NC}: [${align_done}] [${markdup_done}] [${variant_done}] - ${CYAN}STARTED${NC} ($latest_time)"
            else
                echo -e "    ${sample}: [${align_done}] [${markdup_done}] [${variant_done}] - PENDING"
            fi
            
            ((shown++))
        fi
    done
    
    if [[ $shown -eq 0 ]]; then
        echo -e "  ${YELLOW}No samples found${NC}"
    fi
    
    echo ""
}

################################################################################
# MAIN MONITOR LOOP
################################################################################

REFRESH_INTERVAL=10
SHOW_SAMPLES=false
COMPACT_MODE=false

show_usage() {
    cat << EOF
${BOLD}SNP Pipeline Monitor - SLURM Compatible${NC}

Usage: $0 [OPTIONS]

Options:
  -i, --interval SECONDS   Refresh interval (default: 10)
  -s, --samples           Show detailed sample status
  -c, --compact           Compact mode (less output)
  -h, --help              Show this help message

${BOLD}Note:${NC} This monitor reads log files and checkpoints only.

${BOLD}Example:${NC}
  # Basic monitoring
  $0

  # Fast refresh with sample details
  $0 -i 5 -s

  # Compact view
  $0 --compact
EOF
}

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--interval)
            REFRESH_INTERVAL="$2"
            shift 2
            ;;
        -s|--samples)
            SHOW_SAMPLES=true
            shift
            ;;
        -c|--compact)
            COMPACT_MODE=true
            shift
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

while true; do
    print_header
    show_pipeline_status
    show_storage_status
    show_sample_progress
    
    if ! $COMPACT_MODE; then
        show_joint_genotyping
        show_recent_logs
    fi
    
    if $SHOW_SAMPLES; then
        show_sample_details
    fi
    
    echo -e "${BOLD}────────────────────────────────────────────────────────────────────────────${NC}"
    echo -e "Refreshing every ${REFRESH_INTERVAL}s | Press Ctrl+C to exit | Use -h for help"
    
    sleep $REFRESH_INTERVAL
done
