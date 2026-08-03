#!/bin/bash

################################################################################
# SNP Variant Calling Pipeline - PIXY VERSION (WITH INVARIANT SITES)
# Features:
# - Generates VCF files with BOTH variant AND invariant sites for pixy
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
WORK_DIR="/home/dmpachon/jorge/TATIANA/snp_calling/snp_calling_pixy"
OUT_DIR="${WORK_DIR}/results"
LOG_DIR="${WORK_DIR}/logs_parallel"
TMP_DIR="${WORK_DIR}/tmp"
CHECKPOINT_DIR="${WORK_DIR}/checkpoints"

# Create directory structure
mkdir -p "${OUT_DIR}"/{alignment,markdup,metrics,variants/{raw,filtered,all_sites}}
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
# VARIANT CALLING - WITH INVARIANT SITES FOR PIXY
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

    log "  Calling variants for ${sample} (GVCF mode with all sites)..."
    
    check_disk_space

    # CRITICAL FOR PIXY: Use -ERC GVCF mode to emit all sites
    gatk --java-options "-Xmx${MEM_PER_SAMPLE}" HaplotypeCaller \
        -R "${REF_GENOME}" \
        -I "${input_bam}" \
        -O "${output_gvcf}" \
        -ERC GVCF \
        --min-base-quality-score ${MIN_BQ} \
        --tmp-dir "${TMP_DIR}/${sample}" \
        2>&1 | tee "${LOG_DIR}/${sample}_haplotypecaller.log"

    create_checkpoint "${sample}" "variant_calling"

    # CLEANUP: Remove markdup BAM if configured
    if ! $KEEP_MARKDUP_BAM && $CLEANUP_ON_SUCCESS; then
        log "  Removing markdup BAM..."
        rm -f "${input_bam}" "${input_bam}.bai"
    fi
}

################################################################################
# PROCESS SAMPLE
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
# JOINT GENOTYPING - WITH ALL SITES FOR PIXY
################################################################################

joint_genotyping() {
    if checkpoint_exists "cohort" "joint_genotyping" && ! $FORCE_RERUN; then
        info "Joint genotyping already completed, skipping..."
        return 0
    fi

    log "Performing joint genotyping with ALL SITES (for pixy)..."

    local combined_gvcf="${OUT_DIR}/variants/raw/cohort.g.vcf.gz"
    local all_sites_vcf="${OUT_DIR}/variants/all_sites/cohort.allsites.vcf.gz"

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

    # GenotypeGVCFs with ALL SITES
    if [[ ! -f "${all_sites_vcf}" ]]; then
        log "  Genotyping combined GVCF (including invariant sites)..."
        gatk --java-options "-Xmx200G" GenotypeGVCFs \
            -R "${REF_GENOME}" \
            -V "${combined_gvcf}" \
            -O "${all_sites_vcf}" \
            --include-non-variant-sites true \
            2>&1 | tee "${LOG_DIR}/genotype_gvcfs_allsites.log"
        
        # Index the VCF
        bcftools index -t "${all_sites_vcf}"
    fi

    create_checkpoint "cohort" "joint_genotyping"
    
    log "  All-sites VCF created: ${all_sites_vcf}"
    info "  This VCF contains BOTH variant and invariant sites for pixy analysis"
}

################################################################################
# VARIANT FILTERING - CREATES FILTERED VCFS BUT KEEPS ALL-SITES VERSION
################################################################################
filter_variants() {
    if checkpoint_exists "cohort" "filtering" && ! $FORCE_RERUN; then
        info "Variant filtering already completed, skipping..."
        return 0
    fi

    log "Filtering variants..."

    local all_sites_vcf="${OUT_DIR}/variants/all_sites/cohort.allsites.vcf.gz"
    local snp_vcf="${OUT_DIR}/variants/filtered/cohort.snps.vcf.gz"
    local filtered_snp_vcf="${OUT_DIR}/variants/filtered/cohort.snps.filtered.vcf.gz"
    local pass_snp_vcf="${OUT_DIR}/variants/filtered/cohort.snps.pass.vcf.gz"
    local final_variant_vcf="${OUT_DIR}/variants/filtered/cohort.snps.final.vcf.gz"
    
    local invariant_vcf="${OUT_DIR}/variants/filtered/cohort.invariant.vcf.gz"
    local filtered_invariant_vcf="${OUT_DIR}/variants/filtered/cohort.invariant.filtered.vcf.gz"
    local pass_invariant_vcf="${OUT_DIR}/variants/filtered/cohort.invariant.pass.vcf.gz"
    local final_invariant_vcf="${OUT_DIR}/variants/filtered/cohort.invariant.final.vcf.gz"
    
    local merged_allsites_vcf="${OUT_DIR}/variants/filtered/cohort.allsites.filtered.vcf.gz"

    # ========================================================================
    # VARIANT SITES FILTERING
    # ========================================================================
    
    # Select SNPs only
    if [[ ! -f "${snp_vcf}" ]]; then
        log "  Selecting SNPs only..."
        gatk SelectVariants \
            -R "${REF_GENOME}" \
            -V "${all_sites_vcf}" \
            -select-type SNP \
            -O "${snp_vcf}" \
            2>&1 | tee "${LOG_DIR}/select_snps.log"
    fi

    # Apply hard filters to variant sites
    if [[ ! -f "${filtered_snp_vcf}" ]]; then
        log "  Applying hard filters to variant sites..."
        gatk VariantFiltration \
            -R "${REF_GENOME}" \
            -V "${snp_vcf}" \
            -O "${filtered_snp_vcf}" \
            --filter-name "QD_filter" --filter-expression "QD < 2.0" \
            --filter-name "FS_filter" --filter-expression "FS > 60.0" \
            --filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
            --filter-name "SOR_filter" --filter-expression "SOR > 3.0" \
            --filter-name "MQRankSum_filter" --filter-expression "MQRankSum < -12.5" \
            --filter-name "ReadPosRankSum_filter" --filter-expression "ReadPosRankSum < -8.0" \
            2>&1 | tee "${LOG_DIR}/filter_snps.log"
    fi

    # Extract PASS variants
    if [[ ! -f "${pass_snp_vcf}" ]]; then
        log "  Extracting PASS variant sites..."
        bcftools view \
            -f PASS \
            -O z \
            -o "${pass_snp_vcf}" \
            "${filtered_snp_vcf}"
        bcftools index -t "${pass_snp_vcf}"
    fi

    # Apply population filters to variant sites
    if [[ ! -f "${final_variant_vcf}" ]]; then
        log "  Applying population filters to variant sites..."
        vcftools --gzvcf "${pass_snp_vcf}" \
            --remove-indels \
            --hwe 0.01 \
            --maf 0.01 \
            --max-missing 0.8 \
            --min-meanDP 20 \
            --max-meanDP 500 \
            --recode --stdout | bgzip -c > "${final_variant_vcf}"
        bcftools index -t "${final_variant_vcf}"
    fi

    # ========================================================================
    # INVARIANT SITES FILTERING
    # ========================================================================
    
    # Select invariant sites only (sites with no alternate allele)
    if [[ ! -f "${invariant_vcf}" ]]; then
        log "  Selecting invariant sites only..."
        bcftools view \
            -e 'ALT!="."' \
            -O z \
            -o "${invariant_vcf}" \
            "${all_sites_vcf}"
        bcftools index -t "${invariant_vcf}"
    fi

    # Apply hard filters to invariant sites
    if [[ ! -f "${filtered_invariant_vcf}" ]]; then
        log "  Applying hard filters to invariant sites..."
        gatk VariantFiltration \
            -R "${REF_GENOME}" \
            -V "${invariant_vcf}" \
            -O "${filtered_invariant_vcf}" \
            --filter-name "DP_filter" --filter-expression "DP < 20" \
            --filter-name "DP_max_filter" --filter-expression "DP > 500" \
            --filter-name "MQ_filter" --filter-expression "MQ < 40.0" \
            2>&1 | tee "${LOG_DIR}/filter_invariant.log"
    fi

    # Extract PASS invariant sites
    if [[ ! -f "${pass_invariant_vcf}" ]]; then
        log "  Extracting PASS invariant sites..."
        bcftools view \
            -f PASS \
            -O z \
            -o "${pass_invariant_vcf}" \
            "${filtered_invariant_vcf}"
        bcftools index -t "${pass_invariant_vcf}"
    fi

    # Apply population filters to invariant sites
    if [[ ! -f "${final_invariant_vcf}" ]]; then
        log "  Applying population filters to invariant sites..."
        vcftools --gzvcf "${pass_invariant_vcf}" \
            --max-missing 0.8 \
            --min-meanDP 20 \
            --max-meanDP 500 \
            --recode --stdout | bgzip -c > "${final_invariant_vcf}"
        bcftools index -t "${final_invariant_vcf}"
    fi

    # ========================================================================
    # MERGE FILTERED VARIANT AND INVARIANT SITES
    # ========================================================================
    
    if [[ ! -f "${merged_allsites_vcf}" ]]; then
        log "  Merging filtered variant and invariant sites..."
        
        # Concatenate variant and invariant VCFs
        bcftools concat \
            --allow-overlaps \
            -O z \
            -o "${merged_allsites_vcf}.tmp" \
            "${final_variant_vcf}" "${final_invariant_vcf}"
        
        # Sort the merged VCF by position
        bcftools sort \
            -O z \
            -o "${merged_allsites_vcf}" \
            "${merged_allsites_vcf}.tmp"
        
        # Index the final merged file
        bcftools index -t "${merged_allsites_vcf}"
        
        # Clean up temporary file
        rm -f "${merged_allsites_vcf}.tmp"
        
        log "  Created merged all-sites VCF: ${merged_allsites_vcf}"
    fi

    # ========================================================================
    # GENERATE STATISTICS
    # ========================================================================
    
    log "  Generating variant statistics..."
    bcftools stats "${all_sites_vcf}" > "${OUT_DIR}/variants/all_sites/cohort.allsites.stats.txt"
    bcftools stats "${final_variant_vcf}" > "${OUT_DIR}/variants/filtered/cohort.snps.final.stats.txt"
    bcftools stats "${final_invariant_vcf}" > "${OUT_DIR}/variants/filtered/cohort.invariant.final.stats.txt"
    bcftools stats "${merged_allsites_vcf}" > "${OUT_DIR}/variants/filtered/cohort.allsites.filtered.stats.txt"

    # ========================================================================
    # SUMMARY
    # ========================================================================
    
    log "  Filtering complete. Output files:"
    log "    - Original all sites:        ${all_sites_vcf}"
    log "    - Filtered variant sites:    ${final_variant_vcf}"
    log "    - Filtered invariant sites:  ${final_invariant_vcf}"
    log "    - Merged filtered all sites: ${merged_allsites_vcf}"

    create_checkpoint "cohort" "filtering"
}

################################################################################
# PIXY PREPARATION - CREATE POPULATIONS FILE TEMPLATE
################################################################################

create_pixy_populations_template() {
    log "Creating pixy populations file template..."

    local pop_file="${OUT_DIR}/variants/all_sites/populations.txt"
    local all_sites_vcf="${OUT_DIR}/variants/all_sites/cohort.allsites.vcf.gz"

    # Extract sample names from VCF
    bcftools query -l "${all_sites_vcf}" > "${OUT_DIR}/variants/all_sites/sample_list.txt"

    # Create template populations file
    cat > "${pop_file}" << 'EOF'
# Pixy populations file
# Format: sample_name<TAB>population_name
# Edit this file to assign samples to populations
# Example:
# sample1       pop1
# sample2       pop1
# sample3       pop2
# sample4       pop2

EOF

    # Add all samples with placeholder population
    while IFS= read -r sample; do
        echo -e "${sample}\tpopulation_placeholder" >> "${pop_file}"
    done < "${OUT_DIR}/variants/all_sites/sample_list.txt"

    log "  Created populations template: ${pop_file}"
    warn "  IMPORTANT: Edit ${pop_file} to assign samples to correct populations before running pixy!"
}


################################################################################
# PIXY USAGE INSTRUCTIONS
################################################################################

print_pixy_instructions() {
    cat << 'EOF'

================================================================================
                        PIXY ANALYSIS INSTRUCTIONS
================================================================================

Your all-sites VCF file is ready for pixy analysis!

Location: ${OUT_DIR}/variants/all_sites/cohort.allsites.vcf.gz

NEXT STEPS:

1. Edit the populations file:
   ${OUT_DIR}/variants/all_sites/populations.txt
   
   Assign each sample to the correct population.

2. Run pixy to calculate population statistics:

   # Calculate pi (nucleotide diversity)
   pixy --stats pi \
        --vcf ${OUT_DIR}/variants/all_sites/cohort.allsites.vcf.gz \
        --populations ${OUT_DIR}/variants/all_sites/populations.txt \
        --window_size 10000 \
        --n_cores 10 \
        --output_folder ${OUT_DIR}/pixy_results \
        --output_prefix pixy_pi

   # Calculate dxy (between-population diversity)
   pixy --stats dxy \
        --vcf ${OUT_DIR}/variants/all_sites/cohort.allsites.vcf.gz \
        --populations ${OUT_DIR}/variants/all_sites/populations.txt \
        --window_size 10000 \
        --n_cores 10 \
        --output_folder ${OUT_DIR}/pixy_results \
        --output_prefix pixy_dxy

   # Calculate fst
   pixy --stats fst \
        --vcf ${OUT_DIR}/variants/all_sites/cohort.allsites.vcf.gz \
        --populations ${OUT_DIR}/variants/all_sites/populations.txt \
        --window_size 10000 \
        --n_cores 10 \
        --output_folder ${OUT_DIR}/pixy_results \
        --output_prefix pixy_fst

   # Or calculate all statistics at once:
   pixy --stats pi fst dxy \
        --vcf ${OUT_DIR}/variants/all_sites/cohort.allsites.vcf.gz \
        --populations ${OUT_DIR}/variants/all_sites/populations.txt \
        --window_size 10000 \
        --n_cores 10 \
        --output_folder ${OUT_DIR}/pixy_results \
        --output_prefix pixy_allstats

IMPORTANT NOTES:
- The VCF contains BOTH variant AND invariant sites (required by pixy)
- Adjust --window_size based on your genome size and analysis needs
- For chromosome/scaffold-specific analysis, use --chromosomes flag
- For site-specific filtering, consider using bcftools before pixy

For more information, visit: https://github.com/ksamuk/pixy

================================================================================
EOF
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
    echo "=== SNP Calling Pipeline (PIXY VERSION) Started: $(date) ===" > "${MAIN_LOG}"
    
    log "Starting SNP calling pipeline WITH INVARIANT SITES for pixy"
    log "Storage-aware mode: Cleaning intermediate files = ${CLEANUP_ON_SUCCESS}"
    log "Rerun detection: Skip completed samples = ${SKIP_COMPLETED_SAMPLES}"
    log "Main log file: ${MAIN_LOG}"
    log "Output directory: ${OUT_DIR}"
    
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
    
    # Joint genotyping with all sites and filtering
    joint_genotyping
    filter_variants
    
    # Create pixy populations template
    create_pixy_populations_template
    
    log "Pipeline complete!"
    log "All-sites VCF for pixy: ${OUT_DIR}/variants/all_sites/cohort.allsites.vcf.gz"
    log "Filtered SNPs VCF: ${OUT_DIR}/variants/filtered/cohort.snps.pass.vcf.gz"
    
    log "  Filtering complete. Output files:"
    log "    - Original all sites:        ${OUT_DIR}/variants/all_sites/cohort.allsites.vcf.gz"
    log "    - Filtered variant sites:    ${OUT_DIR}/variants/filtered/cohort.snps.final.vcf.gz"
    log "    - Filtered invariant sites:  ${OUT_DIR}/variants/filtered/cohort.invariant.final.vcf.gz"
    log "    - Merged filtered all sites: ${OUT_DIR}/variants/filtered/cohort.allsites.filtered.vcf.gz"
 
    final_usage=$(get_disk_usage)
    info "Final disk usage: ${final_usage}%"
    
    # Print pixy usage instructions
    print_pixy_instructions
}

main "$@"
