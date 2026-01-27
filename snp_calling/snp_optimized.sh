#!/bin/bash

################################################################################
# SNP Variant Calling Pipeline - OPTIMIZED FOR STORAGE
# Features:
# - Automatic cleanup of intermediate files
# - Smart rerun detection (skip completed samples)
# - Progress checkpoint system
# - Compressed intermediate files
################################################################################

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
    exit 1
}

warn() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

################################################################################
# CONFIGURATION
################################################################################

# Input/Output paths
READS_DIR="/home/dmpachon/jorge/TATIANA/qc/qc_results/trimmed_reads"
REF_GENOME="/home/dmpachon/jorge/TATIANA/CallicarpaGenome/car_asm.fa"
WORK_DIR="/home/dmpachon/jorge/TATIANA/snp_calling/"
OUT_DIR="${WORK_DIR}/results"
LOG_DIR="${WORK_DIR}/logs_parallel"
TMP_DIR="${WORK_DIR}/tmp"
CHECKPOINT_DIR="${WORK_DIR}/checkpoints"

# Create directory structure
mkdir -p "${OUT_DIR}"/{alignment,markdup,metrics,variants/{raw,filtered}}
mkdir -p "${LOG_DIR}" "${TMP_DIR}" "${CHECKPOINT_DIR}"

# Computational resources
TOTAL_CORES=80
TOTAL_MEM=600  # GB

# Parallel processing configuration
MAX_PARALLEL_SAMPLES=8
THREADS_PER_SAMPLE=10
MEM_PER_SAMPLE="70G"

# Quality thresholds
MIN_MQ=20
MIN_BQ=20

# STORAGE OPTIMIZATION OPTIONS
KEEP_UNSORTED_BAM=false          # Delete unsorted BAM files
KEEP_SORTED_BAM=false            # Delete sorted BAM after marking duplicates
KEEP_MARKDUP_BAM=true            # Keep duplicate-marked BAM (needed for variant calling)
COMPRESS_INTERMEDIATE=true       # Use additional compression for BAM files
CLEANUP_ON_SUCCESS=true          # Remove intermediate files after successful completion

# CHECKPOINT/RERUN OPTIONS
FORCE_RERUN=false                # Set to true to ignore checkpoints and reprocess everything
SKIP_COMPLETED_SAMPLES=true      # Skip samples that have completed checkpoints

################################################################################
# CHECKPOINT MANAGEMENT
################################################################################

checkpoint_exists() {
    local sample=$1
    local step=$2
    local checkpoint_file="${CHECKPOINT_DIR}/${sample}.${step}.done"
    [[ -f "$checkpoint_file" ]]
}

create_checkpoint() {
    local sample=$1
    local step=$2
    local checkpoint_file="${CHECKPOINT_DIR}/${sample}.${step}.done"
    date > "$checkpoint_file"
    log "  ✓ Checkpoint created: ${sample}.${step}"
}

remove_checkpoint() {
    local sample=$1
    local step=$2
    local checkpoint_file="${CHECKPOINT_DIR}/${sample}.${step}.done"
    rm -f "$checkpoint_file"
}

get_disk_usage() {
    df -h "${WORK_DIR}" | awk 'NR==2 {print $5}' | sed 's/%//'
}

check_disk_space() {
    local usage=$(get_disk_usage)
    if [[ $usage -gt 99 ]]; then
        error "Disk usage is ${usage}%. Please free up space before continuing."
    elif [[ $usage -gt 80 ]]; then
        warn "Disk usage is ${usage}%. Consider freeing up space."
    fi
}

################################################################################
# DEPENDENCY CHECK
################################################################################

check_dependencies() {
    log "Checking dependencies..."
    local deps=("samtools" "gatk" "bcftools" "parallel")
    for dep in "${deps[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            error "$dep is not installed or not in PATH"
        fi
    done
    log "All dependencies found"
}

################################################################################
# REFERENCE GENOME PREPARATION
################################################################################

prepare_reference() {
    log "Preparing reference genome..."

    # Index reference with BWA
    if [[ ! -f "${REF_GENOME}.bwt" ]]; then
        log "Creating BWA index..."
        bwa index "${REF_GENOME}" 2>&1 | tee "${LOG_DIR}/bwa_index.log"
    else
        log "BWA index already exists"
    fi

    # Create FASTA index
    if [[ ! -f "${REF_GENOME}.fai" ]]; then
        log "Creating FASTA index..."
        samtools faidx "${REF_GENOME}"
    else
        log "FASTA index already exists"
    fi

    # Create sequence dictionary
    local dict="${REF_GENOME%.fa}.dict"
    if [[ ! -f "${dict}" ]]; then
        log "Creating sequence dictionary..."
        gatk CreateSequenceDictionary \
            -R "${REF_GENOME}" \
            -O "${dict}" 2>&1 | tee "${LOG_DIR}/dict.log"
    else
        log "Sequence dictionary already exists"
    fi
}

################################################################################
# READ ALIGNMENT
################################################################################

align_reads() {
    local sample=$1
    local r1=$2
    local r2=$3

    # Check if already completed
    if $SKIP_COMPLETED_SAMPLES && checkpoint_exists "${sample}" "alignment" && ! $FORCE_RERUN; then
        info "  Alignment already completed for ${sample}, skipping..."
        return 0
    fi

    log "  Aligning reads for ${sample}..."

    local sorted_bam="${OUT_DIR}/alignment/${sample}.sorted.bam"

    check_disk_space

    # Alignment with BWA-MEM, piped directly to sorting
    bwa mem \
        -t ${THREADS_PER_SAMPLE} \
        -M \
        -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
        "${REF_GENOME}" \
        "${r1}" "${r2}" \
        2> "${LOG_DIR}/${sample}_bwa.log" \
        | samtools sort \
            -@ ${THREADS_PER_SAMPLE} \
            -m 4G \
            -l 9 \
            -o "${sorted_bam}" \
            - 2>&1 | tee -a "${LOG_DIR}/${sample}_sort.log"

    # Index BAM
    samtools index -@ ${THREADS_PER_SAMPLE} "${sorted_bam}"

    # Generate alignment statistics (lightweight)
    samtools flagstat \
        -@ ${THREADS_PER_SAMPLE} \
        "${sorted_bam}" > "${OUT_DIR}/metrics/${sample}_flagstat.txt"

    create_checkpoint "${sample}" "alignment"
}

################################################################################
# MARK DUPLICATES
################################################################################

mark_duplicates() {
    local sample=$1
    
    # Check if already completed
    if $SKIP_COMPLETED_SAMPLES && checkpoint_exists "${sample}" "markdup" && ! $FORCE_RERUN; then
        info "  Duplicate marking already completed for ${sample}, skipping..."
        return 0
    fi

    local input_bam="${OUT_DIR}/alignment/${sample}.sorted.bam"
    local output_bam="${OUT_DIR}/markdup/${sample}.markdup.bam"
    local metrics="${OUT_DIR}/metrics/${sample}_duplicate_metrics.txt"

    log "  Marking duplicates for ${sample}..."
    
    check_disk_space

    gatk --java-options "-Xmx${MEM_PER_SAMPLE}" MarkDuplicates \
        -I "${input_bam}" \
        -O "${output_bam}" \
        -M "${metrics}" \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT \
        --COMPRESSION_LEVEL 9 \
        --TMP_DIR "${TMP_DIR}/${sample}" \
        2>&1 | tee "${LOG_DIR}/${sample}_markdup.log"

    create_checkpoint "${sample}" "markdup"

    # CLEANUP: Remove sorted BAM if configured
    if ! $KEEP_SORTED_BAM; then
        log "  Removing intermediate sorted BAM..."
        rm -f "${input_bam}" "${input_bam}.bai"
    fi
}

################################################################################
# VARIANT CALLING
################################################################################

call_variants_haplotypecaller() {
    local sample=$1
    
    # Check if already completed
    if $SKIP_COMPLETED_SAMPLES && checkpoint_exists "${sample}" "variant_calling" && ! $FORCE_RERUN; then
        info "  Variant calling already completed for ${sample}, skipping..."
        return 0
    fi

    local input_bam="${OUT_DIR}/markdup/${sample}.markdup.bam"
    local output_gvcf="${OUT_DIR}/variants/raw/${sample}.g.vcf.gz"

    log "  Calling variants for ${sample}..."
    
    check_disk_space

    gatk --java-options "-Xmx${MEM_PER_SAMPLE}" HaplotypeCaller \
        -R "${REF_GENOME}" \
        -I "${input_bam}" \
        -O "${output_gvcf}" \
        -ERC GVCF \
        --native-pair-hmm-threads ${THREADS_PER_SAMPLE} \
        --min-base-quality-score ${MIN_BQ} \
        2>&1 | tee "${LOG_DIR}/${sample}_haplotypecaller.log"

    create_checkpoint "${sample}" "variant_calling"

    # CLEANUP: Remove markdup BAM if configured (only after successful GVCF creation)
    if ! $KEEP_MARKDUP_BAM && [[ -f "${output_gvcf}" ]]; then
        warn "  Removing markdup BAM (can be regenerated from sorted BAM if needed)..."
        rm -f "${input_bam}" "${input_bam}.bai"
    fi
}

################################################################################
# PROCESS SINGLE SAMPLE
################################################################################

process_sample() {
    local sample=$1
    local r1=$2
    local r2=$3
    
    # Check if sample is fully completed
    if $SKIP_COMPLETED_SAMPLES && \
       checkpoint_exists "${sample}" "alignment" && \
       checkpoint_exists "${sample}" "markdup" && \
       checkpoint_exists "${sample}" "variant_calling" && \
       ! $FORCE_RERUN; then
        info "Sample ${sample} already fully processed, skipping all steps..."
        return 0
    fi
    
    mkdir -p "${TMP_DIR}/${sample}"
    
    log "Starting processing for ${sample}"
    
    align_reads "${sample}" "${r1}" "${r2}"
    mark_duplicates "${sample}"
    call_variants_haplotypecaller "${sample}"
    
    # Cleanup sample tmp directory
    rm -rf "${TMP_DIR}/${sample}"
    
    log "Completed processing for ${sample}"
    
    # Show current disk usage
    local usage=$(get_disk_usage)
    info "Current disk usage: ${usage}%"
}

# Export functions and variables for GNU parallel
export -f process_sample align_reads mark_duplicates call_variants_haplotypecaller
export -f log warn error info checkpoint_exists create_checkpoint get_disk_usage check_disk_space
export OUT_DIR LOG_DIR TMP_DIR CHECKPOINT_DIR REF_GENOME 
export THREADS_PER_SAMPLE MEM_PER_SAMPLE MIN_BQ
export KEEP_SORTED_BAM KEEP_MARKDUP_BAM SKIP_COMPLETED_SAMPLES FORCE_RERUN
export GREEN RED YELLOW BLUE NC WORK_DIR

################################################################################
# JOINT GENOTYPING
################################################################################

joint_genotyping() {
    if checkpoint_exists "cohort" "joint_genotyping" && ! $FORCE_RERUN; then
        info "Joint genotyping already completed, skipping..."
        return 0
    fi

    log "Performing joint genotyping..."

    local combined_gvcf="${OUT_DIR}/variants/raw/cohort.g.vcf.gz"
    local raw_vcf="${OUT_DIR}/variants/raw/cohort.raw.vcf.gz"

    # Collect GVCF files
    local gvcf_files=()
    for gvcf in "${OUT_DIR}"/variants/raw/*.g.vcf.gz; do
        if [[ -f "$gvcf" ]] && [[ "$gvcf" != *"cohort"* ]]; then
            gvcf_files+=("-V" "$gvcf")
        fi
    done

    if [[ ${#gvcf_files[@]} -eq 0 ]]; then
        error "No GVCF files found for joint genotyping"
    fi

    # Combine GVCFs
    if [[ ! -f "${combined_gvcf}" ]]; then
        log "  Combining ${#gvcf_files[@]} GVCFs..."
        gatk --java-options "-Xmx100G" CombineGVCFs \
            -R "${REF_GENOME}" \
            "${gvcf_files[@]}" \
            -O "${combined_gvcf}" \
            2>&1 | tee "${LOG_DIR}/combine_gvcfs.log"
    fi

    # GenotypeGVCFs
    if [[ ! -f "${raw_vcf}" ]]; then
        log "  Genotyping combined GVCF..."
        gatk --java-options "-Xmx200G" GenotypeGVCFs \
            -R "${REF_GENOME}" \
            -V "${combined_gvcf}" \
            -O "${raw_vcf}" \
            2>&1 | tee "${LOG_DIR}/genotype_gvcfs.log"
    fi

    create_checkpoint "cohort" "joint_genotyping"
}

################################################################################
# VARIANT FILTERING
################################################################################

filter_variants() {
    if checkpoint_exists "cohort" "filtering" && ! $FORCE_RERUN; then
        info "Variant filtering already completed, skipping..."
        return 0
    fi

    log "Filtering variants..."

    local raw_vcf="${OUT_DIR}/variants/raw/cohort.raw.vcf.gz"
    local snp_vcf="${OUT_DIR}/variants/filtered/cohort.snps.vcf.gz"
    local filtered_vcf="${OUT_DIR}/variants/filtered/cohort.snps.filtered.vcf.gz"
    local pass_vcf="${OUT_DIR}/variants/filtered/cohort.snps.pass.vcf.gz"

    # Select SNPs
    if [[ ! -f "${snp_vcf}" ]]; then
        log "  Selecting SNPs..."
        gatk SelectVariants \
            -R "${REF_GENOME}" \
            -V "${raw_vcf}" \
            -select-type SNP \
            -O "${snp_vcf}" \
            2>&1 | tee "${LOG_DIR}/select_snps.log"
    fi

    # Apply filters
    if [[ ! -f "${filtered_vcf}" ]]; then
        log "  Applying hard filters..."
        gatk VariantFiltration \
            -R "${REF_GENOME}" \
            -V "${snp_vcf}" \
            -O "${filtered_vcf}" \
            --filter-name "QD_filter" --filter-expression "QD < 2.0" \
            --filter-name "FS_filter" --filter-expression "FS > 60.0" \
            --filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
            --filter-name "SOR_filter" --filter-expression "SOR > 3.0" \
            --filter-name "MQRankSum_filter" --filter-expression "MQRankSum < -12.5" \
            --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0" \
            2>&1 | tee "${LOG_DIR}/filter_snps.log"
    fi

    # Extract PASS variants
    if [[ ! -f "${pass_vcf}" ]]; then
        log "  Extracting PASS variants..."
        bcftools view \
            -f PASS \
            -O z \
            -o "${pass_vcf}" \
            "${filtered_vcf}"
        bcftools index -t "${pass_vcf}"
    fi

    # Generate statistics
    log "  Generating variant statistics..."
    bcftools stats "${pass_vcf}" > "${OUT_DIR}/variants/filtered/cohort.snps.pass.stats.txt"

    create_checkpoint "cohort" "filtering"
}

################################################################################
# MAIN PIPELINE
################################################################################

# Setup main log file for remote monitoring
MAIN_LOG="${LOG_DIR}/pipeline_main.log"

# Override log function to also write to main log
log() {
    local msg="${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
    echo -e "$msg"
    echo -e "$msg" | sed 's/\x1b\[[0-9;]*m//g' >> "${MAIN_LOG}"
}

error() {
    local msg="${RED}[ERROR]${NC} $1"
    echo -e "$msg" >&2
    echo -e "$msg" | sed 's/\x1b\[[0-9;]*m//g' >> "${MAIN_LOG}"
    exit 1
}

warn() {
    local msg="${YELLOW}[WARNING]${NC} $1"
    echo -e "$msg"
    echo -e "$msg" | sed 's/\x1b\[[0-9;]*m//g' >> "${MAIN_LOG}"
}

info() {
    local msg="${BLUE}[INFO]${NC} $1"
    echo -e "$msg"
    echo -e "$msg" | sed 's/\x1b\[[0-9;]*m//g' >> "${MAIN_LOG}"
}

main() {
    # Initialize main log
    echo "=== SNP Calling Pipeline Started: $(date) ===" > "${MAIN_LOG}"
    
    log "Starting OPTIMIZED SNP calling pipeline"
    log "Storage-aware mode: Cleaning intermediate files = ${CLEANUP_ON_SUCCESS}"
    log "Rerun detection: Skip completed samples = ${SKIP_COMPLETED_SAMPLES}"
    log "Main log file: ${MAIN_LOG}"
    
    check_dependencies
    check_disk_space
    prepare_reference
    
    # Discover samples
    log "Discovering sample pairs..."
    declare -A samples
    
    for r1 in "${READS_DIR}"/*_R1*.fastq.gz "${READS_DIR}"/*_1.fastq.gz; do
        if [[ -f "$r1" ]]; then
            sample=$(basename "$r1" | sed -E 's/_R?1.*\.fastq\.gz$//')
            r2="${r1/_R1/_R2}"
            r2="${r2/_1./_2.}"
            
            if [[ -f "$r2" ]]; then
                samples["$sample"]="$r1:$r2"
                
                # Check checkpoint status
                if checkpoint_exists "${sample}" "variant_calling" && $SKIP_COMPLETED_SAMPLES; then
                    info "  ${sample} (COMPLETED ✓)"
                else
                    log "  ${sample} (TO PROCESS)"
                fi
            fi
        fi
    done
    
    if [[ ${#samples[@]} -eq 0 ]]; then
        error "No sample pairs found in ${READS_DIR}"
    fi
    
    log "Found ${#samples[@]} samples total"
    
    # Create job list
    job_file="${TMP_DIR}/jobs.txt"
    > "$job_file"
    
    for sample in "${!samples[@]}"; do
        IFS=':' read -r r1 r2 <<< "${samples[$sample]}"
        echo "process_sample \"${sample}\" \"${r1}\" \"${r2}\"" >> "$job_file"
    done
    
    # Process in parallel
    log "Processing samples (max ${MAX_PARALLEL_SAMPLES} parallel)..."
    parallel --jobs ${MAX_PARALLEL_SAMPLES} --bar < "$job_file"
    
    log "All samples processed!"
    
    # Joint genotyping and filtering
    joint_genotyping
    filter_variants
    
    log "Pipeline complete!"
    log "Final VCF: ${OUT_DIR}/variants/filtered/cohort.snps.pass.vcf.gz"
    
    final_usage=$(get_disk_usage)
    info "Final disk usage: ${final_usage}%"
}

main "$@"
