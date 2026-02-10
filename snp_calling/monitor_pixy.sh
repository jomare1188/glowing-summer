#!/bin/bash

################################################################################
# SNP Pipeline Monitor - PIXY VERSION - SLURM Compatible
# Monitors pipeline progress via log files and checkpoint files only
# Updated for pipeline with invariant sites
# Can run on different machine from the main pipeline
################################################################################

# Configuration (match your pipeline settings)
WORK_DIR="/home/dmpachon/jorge/TATIANA/snp_calling/snp_calling_pixy"
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
    echo -e "${BOLD}${CYAN}║          SNP PIPELINE MONITOR - PIXY VERSION (with invariants)            ║${NC}"
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
    echo -e "  ${MAGENTA}Mode: PIXY-compatible (includes invariant sites)${NC}"
    
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
    local all_sites_size=$(get_dir_size "${OUT_DIR}/variants/all_sites")
    local filtered_size=$(get_dir_size "${OUT_DIR}/variants/filtered")
    local logs_size=$(get_dir_size "${LOG_DIR}")
    local total_size=$((alignment_size + markdup_size + variants_size + logs_size))
    
    echo ""
    echo -e "  ${BOLD}Storage breakdown:${NC}"
    echo -e "    Alignment BAMs:     $(human_readable_size $alignment_size)"
    echo -e "    Duplicate-marked:   $(human_readable_size $markdup_size)"
    echo -e "    Variants (total):   $(human_readable_size $variants_size)"
    echo -e "      ├─ All sites:     $(human_readable_size $all_sites_size) ${MAGENTA}(for pixy)${NC}"
    echo -e "      └─ Filtered SNPs: $(human_readable_size $filtered_size)"
    echo -e "    Logs:               $(human_readable_size $logs_size)"
    echo -e "    ${BOLD}Total pipeline:     $(human_readable_size $total_size)${NC}"
    
    if [[ $all_sites_size -gt 0 ]] && [[ $filtered_size -gt 0 ]]; then
        local size_ratio=$(awk "BEGIN {printf \"%.1f\", $all_sites_size / $filtered_size}")
        echo -e "    ${CYAN}All-sites VCF is ${size_ratio}x larger than filtered VCF${NC}"
    fi
    
    echo ""
}

show_sample_progress() {
    echo -e "${BOLD}${YELLOW}[3] SAMPLE PROGRESS${NC}"
    echo -e "${BOLD}────────────────────────────────────────────────────────────────────────────${NC}"
    
    local total_samples=0
    shopt -s nullglob
    for r1 in "${READS_DIR}"/*_R1*.fastq.gz "${READS_DIR}"/*_1.fastq.gz; do
        if [[ -f "$r1" ]]; then
            ((total_samples++))
        fi
    done
    shopt -u nullglob
    
    local aligned=0
    local markdup=0
    local variant_called=0
    
    if [[ -d "${CHECKPOINT_DIR}" ]]; then
        aligned=$(ls "${CHECKPOINT_DIR}"/*.alignment.done 2>/dev/null | wc -l)
        markdup=$(ls "${CHECKPOINT_DIR}"/*.markdup.done 2>/dev/null | wc -l)
        variant_called=$(ls "${CHECKPOINT_DIR}"/*.variant_calling.done 2>/dev/null | wc -l)
    fi
    
    local align_pct=0
    local markdup_pct=0
    local variant_pct=0
    
    if [[ $total_samples -gt 0 ]]; then
        align_pct=$(awk "BEGIN {printf \"%.1f\", ($aligned / $total_samples) * 100}")
        markdup_pct=$(awk "BEGIN {printf \"%.1f\", ($markdup / $total_samples) * 100}")
        variant_pct=$(awk "BEGIN {printf \"%.1f\", ($variant_called / $total_samples) * 100}")
    fi
    
    echo -e "  ${BOLD}Total samples: ${total_samples}${NC}"
    echo ""
    echo -e "  Alignment:        ${aligned}/${total_samples} (${align_pct}%)"
    echo -e "  Mark duplicates:  ${markdup}/${total_samples} (${markdup_pct}%)"
    echo -e "  Variant calling:  ${variant_called}/${total_samples} (${variant_pct}%)"
    
    local gvcf_count=$(count_files "${OUT_DIR}/variants/raw" "*.g.vcf.gz")
    echo -e "  GVCF files:       ${gvcf_count}"
    
    echo ""
    
    if [[ $total_samples -gt 0 ]]; then
        local completed=$variant_called
        local in_progress=$((markdup - variant_called))
        local pending=$((total_samples - markdup))
        
        echo -e "  ${GREEN}Completed:${NC}    ${completed}"
        echo -e "  ${YELLOW}In progress:${NC}  ${in_progress}"
        echo -e "  ${BLUE}Pending:${NC}      ${pending}"
    fi
    
    echo ""
}

show_joint_genotyping() {
    echo -e "${BOLD}${BLUE}[4] JOINT GENOTYPING & FILTERING (PIXY MODE)${NC}"
    echo -e "${BOLD}────────────────────────────────────────────────────────────────────────────${NC}"

    local joint_done=false
    local filter_done=false

    [[ -f "${CHECKPOINT_DIR}/cohort.joint_genotyping.done" ]] && joint_done=true
    [[ -f "${CHECKPOINT_DIR}/cohort.filtering.done" ]] && filter_done=true

    # Files specific to pixy pipeline
    local combined_gvcf="${OUT_DIR}/variants/raw/cohort.g.vcf.gz"
    local all_sites_vcf="${OUT_DIR}/variants/all_sites/cohort.allsites.vcf.gz"
    local filtered_variant_vcf="${OUT_DIR}/variants/filtered/cohort.snps.final.vcf.gz"
    local filtered_invariant_vcf="${OUT_DIR}/variants/filtered/cohort.invariant.final.vcf.gz"
    local merged_allsites_vcf="${OUT_DIR}/variants/filtered/cohort.allsites.filtered.vcf.gz"
    local populations_file="${OUT_DIR}/variants/all_sites/populations.txt"

    echo -e "  ${BOLD}VCF Files:${NC}"
    echo ""

    if [[ -f "$combined_gvcf" ]]; then
        local size=$(stat -c%s "$combined_gvcf" 2>/dev/null || stat -f%z "$combined_gvcf" 2>/dev/null || echo "0")
        local age=$(get_file_age "$combined_gvcf")
        echo -e "  ${GREEN}✓${NC} Combined GVCF:         $(human_readable_size $size) ($age)"
    else
        echo -e "  ${YELLOW}○${NC} Combined GVCF:         Not created yet"
    fi

    if [[ -f "$all_sites_vcf" ]]; then
        local size=$(stat -c%s "$all_sites_vcf" 2>/dev/null || stat -f%z "$all_sites_vcf" 2>/dev/null || echo "0")
        local age=$(get_file_age "$all_sites_vcf")

        local total_sites="unknown"
        local variant_sites="unknown"
        local invariant_sites="unknown"

        if command -v bcftools &> /dev/null; then
            total_sites=$(bcftools view -H "$all_sites_vcf" 2>/dev/null | wc -l | tr -d ' ')
            variant_sites=$(bcftools view -H -v snps "$all_sites_vcf" 2>/dev/null | wc -l | tr -d ' ')
            if [[ "$total_sites" =~ ^[0-9]+$ ]] && [[ "$variant_sites" =~ ^[0-9]+$ ]]; then
                invariant_sites=$((total_sites - variant_sites))
            fi
        fi

        echo -e "  ${GREEN}✓${NC} ${MAGENTA}All-sites VCF (raw):${NC}   $(human_readable_size $size) ($age)"
        if [[ "$total_sites" != "unknown" ]]; then
            local inv_pct=$(awk "BEGIN {printf \"%.1f\", ($invariant_sites / $total_sites) * 100}")
            echo -e "      Total sites:       ${total_sites}"
            echo -e "      Variant sites:     ${variant_sites}"
            echo -e "      Invariant sites:   ${invariant_sites} (${inv_pct}%)"
        fi
    else
        echo -e "  ${YELLOW}○${NC} ${MAGENTA}All-sites VCF (raw):${NC}   Not created yet"
    fi

    echo ""
    echo -e "  ${BOLD}Filtered VCF Files:${NC}"
    echo ""

    if [[ -f "$filtered_variant_vcf" ]]; then
        local size=$(stat -c%s "$filtered_variant_vcf" 2>/dev/null || stat -f%z "$filtered_variant_vcf" 2>/dev/null || echo "0")
        local age=$(get_file_age "$filtered_variant_vcf")

        local variant_count="unknown"
        if command -v bcftools &> /dev/null; then
            variant_count=$(bcftools view -H "$filtered_variant_vcf" 2>/dev/null | wc -l | tr -d ' ')
        fi

        echo -e "  ${GREEN}✓${NC} Filtered variants:     $(human_readable_size $size) (${variant_count} SNPs, $age)"
    else
        echo -e "  ${YELLOW}○${NC} Filtered variants:     Not created yet"
    fi

    if [[ -f "$filtered_invariant_vcf" ]]; then
        local size=$(stat -c%s "$filtered_invariant_vcf" 2>/dev/null || stat -f%z "$filtered_invariant_vcf" 2>/dev/null || echo "0")
        local age=$(get_file_age "$filtered_invariant_vcf")

        local invariant_count="unknown"
        if command -v bcftools &> /dev/null; then
            invariant_count=$(bcftools view -H "$filtered_invariant_vcf" 2>/dev/null | wc -l | tr -d ' ')
        fi

        echo -e "  ${GREEN}✓${NC} Filtered invariants:   $(human_readable_size $size) (${invariant_count} sites, $age)"
    else
        echo -e "  ${YELLOW}○${NC} Filtered invariants:   Not created yet"
    fi

    if [[ -f "$merged_allsites_vcf" ]]; then
        local size=$(stat -c%s "$merged_allsites_vcf" 2>/dev/null || stat -f%z "$merged_allsites_vcf" 2>/dev/null || echo "0")
        local age=$(get_file_age "$merged_allsites_vcf")

        local total_count="unknown"
        local var_count="unknown"
        local inv_count="unknown"
        
        if command -v bcftools &> /dev/null; then
            total_count=$(bcftools view -H "$merged_allsites_vcf" 2>/dev/null | wc -l | tr -d ' ')
            var_count=$(bcftools view -H -v snps "$merged_allsites_vcf" 2>/dev/null | wc -l | tr -d ' ')
            if [[ "$total_count" =~ ^[0-9]+$ ]] && [[ "$var_count" =~ ^[0-9]+$ ]]; then
                inv_count=$((total_count - var_count))
            fi
        fi

        echo -e "  ${GREEN}✓${NC} ${CYAN}Merged all-sites (filtered):${NC} $(human_readable_size $size) ($age)"
        if [[ "$total_count" != "unknown" ]]; then
            echo -e "      Total sites:       ${total_count}"
            echo -e "      Variant sites:     ${var_count}"
            echo -e "      Invariant sites:   ${inv_count}"
            echo -e "      ${GREEN}Ready for pixy!${NC}"
        fi
    else
        echo -e "  ${YELLOW}○${NC} ${CYAN}Merged all-sites (filtered):${NC} Not created yet ${YELLOW}(needed for pixy!)${NC}"
    fi

    echo ""

    # Check populations file
    if [[ -f "$populations_file" ]]; then
        local age=$(get_file_age "$populations_file")
        if grep -q "population_placeholder" "$populations_file" 2>/dev/null; then
            echo -e "  ${YELLOW}⚠${NC} Populations file: Present but ${YELLOW}NOT EDITED${NC} ($age)"
            echo -e "      ${YELLOW}Edit this file before running pixy!${NC}"
        else
            local pop_count=$(tail -n +2 "$populations_file" 2>/dev/null | awk '{print $2}' | sort -u | wc -l)
            echo -e "  ${GREEN}✓${NC} Populations file: Configured with ${pop_count} populations ($age)"
        fi
    else
        echo -e "  ${YELLOW}○${NC} Populations file: Not created yet"
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

    # Pixy readiness check
    echo ""
    if [[ -f "$merged_allsites_vcf" ]] && [[ -f "$populations_file" ]] && ! grep -q "population_placeholder" "$populations_file" 2>/dev/null; then
        echo -e "  ${GREEN}✓✓✓ Ready for pixy analysis!${NC}"
        echo -e "      Input VCF: ${merged_allsites_vcf}"
        echo -e "      Populations: ${populations_file}"
    elif [[ -f "$merged_allsites_vcf" ]]; then
        echo -e "  ${YELLOW}⚠ Almost ready: Edit populations.txt before running pixy${NC}"
    elif [[ -f "$all_sites_vcf" ]]; then
        echo -e "  ${BLUE}○ Filtering in progress or pending${NC}"
    else
        echo -e "  ${BLUE}○ Waiting for all-sites VCF generation${NC}"
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

show_pixy_status() {
    echo -e "${BOLD}${MAGENTA}[7] PIXY ANALYSIS STATUS${NC}"
    echo -e "${BOLD}────────────────────────────────────────────────────────────────────────────${NC}"
    
    local pixy_out="${OUT_DIR}/pixy_results"
    
    if [[ ! -d "$pixy_out" ]]; then
        echo -e "  ${BLUE}○ No pixy results directory found${NC}"
        echo -e "  Run pixy after pipeline completes and populations.txt is edited"
    else
        local pi_files=$(count_files "$pixy_out" "*_pi.txt")
        local dxy_files=$(count_files "$pixy_out" "*_dxy.txt")
        local fst_files=$(count_files "$pixy_out" "*_fst.txt")
        
        echo -e "  ${GREEN}✓${NC} Pixy results directory found"
        echo ""
        echo -e "  Results found:"
        echo -e "    Pi (diversity):       ${pi_files} file(s)"
        echo -e "    Dxy (divergence):     ${dxy_files} file(s)"
        echo -e "    Fst (differentiation): ${fst_files} file(s)"
        
        if [[ $pi_files -gt 0 ]] || [[ $dxy_files -gt 0 ]] || [[ $fst_files -gt 0 ]]; then
            local latest_result=$(ls -t "$pixy_out"/*.txt 2>/dev/null | head -1)
            if [[ -n "$latest_result" ]]; then
                echo ""
                echo -e "  Latest result: $(basename "$latest_result") ($(get_file_age "$latest_result"))"
            fi
        fi
    fi
    
    echo ""
}

################################################################################
# MAIN MONITOR LOOP
################################################################################

REFRESH_INTERVAL=10
SHOW_SAMPLES=true
SHOW_PIXY=true
COMPACT_MODE=false

show_usage() {
    cat << EOF
${BOLD}SNP Pipeline Monitor - PIXY VERSION${NC}

Usage: $0 [OPTIONS]

Options:
  -i, --interval SECONDS   Refresh interval (default: 10)
  -s, --samples           Show detailed sample status
  -p, --pixy              Show pixy analysis status
  -c, --compact           Compact mode (less output)
  -h, --help              Show this help message

${BOLD}Note:${NC} This monitor reads log files and checkpoints only.
${BOLD}Pipeline mode:${NC} PIXY-compatible (includes invariant sites)

${BOLD}Example:${NC}
  # Basic monitoring
  $0

  # Fast refresh with sample and pixy details
  $0 -i 5 -s -p

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
        -p|--pixy)
            SHOW_PIXY=true
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
    
    if $SHOW_PIXY; then
        show_pixy_status
    fi
    
    echo -e "${BOLD}────────────────────────────────────────────────────────────────────────────${NC}"
    echo -e "Refreshing every ${REFRESH_INTERVAL}s | Press Ctrl+C to exit | Use -h for help"
    
    sleep $REFRESH_INTERVAL
done
