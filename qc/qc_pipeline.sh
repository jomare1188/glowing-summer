#!/bin/bash
# Illumina Paired-End WGS Quality Control Pipeline
# Usage: ./qc_pipeline.sh <input_dir> <output_dir> <threads>

set -euo pipefail

# Parameters
INPUT_DIR=${1:-.}
OUTPUT_DIR=${2:-qc_results}
THREADS=${3:-8}

# Create directory structure
mkdir -p ${OUTPUT_DIR}/{raw_fastqc,trimmed_fastqc,trimmed_reads,logs}

echo "Starting QC pipeline..."
echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"

# Step 1: Initial FastQC on raw reads
echo "[$(date)] Running FastQC on raw reads..."
fastqc -t ${THREADS} \
    -o ${OUTPUT_DIR}/raw_fastqc \
    ${INPUT_DIR}/*fastq.gz \
    2>&1 | tee ${OUTPUT_DIR}/logs/raw_fastqc.log

## Step 2: MultiQC report for raw reads
echo "[$(date)] Generating MultiQC report for raw reads..."
multiqc ${OUTPUT_DIR}/raw_fastqc \
    -o ${OUTPUT_DIR}/raw_fastqc \
    -n raw_multiqc_report \
    2>&1 | tee ${OUTPUT_DIR}/logs/raw_multiqc.log

# Step 3: Trimming and filtering with fastp (parallelized)
echo "[$(date)] Running fastp for trimming and filtering..."

# Calculate threads per job
PARALLEL_JOBS=10
THREADS_PER_JOB=$((THREADS / PARALLEL_JOBS))

# Ensure at least 1 thread per job
if [ ${THREADS_PER_JOB} -lt 1 ]; then
    THREADS_PER_JOB=1
fi

echo "Running ${PARALLEL_JOBS} parallel fastp jobs with ${THREADS_PER_JOB} threads each"

# Export variables needed by parallel function
export OUTPUT_DIR
export THREADS_PER_JOB

# Define the processing function
process_fastp() {
    R1=$1
    SAMPLE=$(basename ${R1} | sed 's/_R1.*//')
    R2=${R1/_R1/_R2}
    
    echo "[$(date)] Processing sample: ${SAMPLE}"
    
    fastp \
        -i ${R1} \
        -I ${R2} \
        -o ${OUTPUT_DIR}/trimmed_reads/${SAMPLE}_R1_trimmed.fastq.gz \
        -O ${OUTPUT_DIR}/trimmed_reads/${SAMPLE}_R2_trimmed.fastq.gz \
        --thread ${THREADS_PER_JOB} \
	-f 14 \
	-F 14 \
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
        --html ${OUTPUT_DIR}/trimmed_reads/${SAMPLE}_fastp.html \
        --json ${OUTPUT_DIR}/trimmed_reads/${SAMPLE}_fastp.json \
        2>&1 | tee ${OUTPUT_DIR}/logs/${SAMPLE}_fastp.log
}

# Export the function so parallel can use it
export -f process_fastp

# Run parallel processing
ls ${INPUT_DIR}/*_R1*.fastq.gz | \
    parallel --env process_fastp --env OUTPUT_DIR --env THREADS_PER_JOB \
        -j ${PARALLEL_JOBS} --verbose process_fastp {}

# Step 4: FastQC on trimmed reads
echo "[$(date)] Running FastQC on trimmed reads..."
fastqc -t ${THREADS} \
    -o ${OUTPUT_DIR}/trimmed_fastqc \
    ${OUTPUT_DIR}/trimmed_reads/*_trimmed.fastq.gz \
    2>&1 | tee ${OUTPUT_DIR}/logs/trimmed_fastqc.log

# Step 5: MultiQC report for trimmed reads
echo "[$(date)] Generating MultiQC report for trimmed reads..."
multiqc ${OUTPUT_DIR}/trimmed_fastqc \
    ${OUTPUT_DIR}/trimmed_reads \
    -o ${OUTPUT_DIR}/trimmed_fastqc \
    -n trimmed_multiqc_report \
    2>&1 | tee ${OUTPUT_DIR}/logs/trimmed_multiqc.log

# Step 6: Generate summary statistics
echo "[$(date)] Generating summary statistics..."
{
    echo "QC Pipeline Summary"
    echo "==================="
    echo "Date: $(date)"
    echo "Raw reads processed: $(ls ${INPUT_DIR}/*_R1*.fastq.gz | wc -l)"
    echo "Trimmed reads generated: $(ls ${OUTPUT_DIR}/trimmed_reads/*_R1_trimmed.fastq.gz | wc -l)"
    echo ""
    echo "Check the following reports:"
    echo "- Raw reads: ${OUTPUT_DIR}/raw_fastqc/raw_multiqc_report.html"
    echo "- Trimmed reads: ${OUTPUT_DIR}/trimmed_fastqc/trimmed_multiqc_report.html"
    echo "- Individual fastp reports in: ${OUTPUT_DIR}/trimmed_reads/"
} | tee ${OUTPUT_DIR}/qc_summary.txt

echo "[$(date)] QC pipeline completed successfully!"
echo "Results are in: ${OUTPUT_DIR}"
