#!/bin/sh
#SBATCH --job-name=filter_vcf
#SBATCH --partition=long
#SBATCH --ntasks-per-node=10
#SBATCH --mem=30gb
#SBATCH --error=%j_snp_filter.err
#SBATCH --output=%j_snp_filter.out


source ~/.bashrc
conda activate SNP_call

IN="../cohort.snps.final.vcf.gz"
LD_FILTERED="ld_pruned.snps.vcf.gz"
BIG="ld_pruned_50K.snps.vcf.gz"
SMALL="ld_pruned_5K.snps.vcf.gz"

filter_ld () {
# LD parameters
local WINDOW_SIZE=50
local LD_THRESHOLD=0.2
local STEP_SIZE=10

echo "Starting LD pruning..."

# Step 1: Convert VCF to PLINK format
echo "Converting VCF to PLINK..."
plink --vcf $1 \
      --make-bed \
      --out temp_data \
      --allow-extra-chr \
      --double-id

# Step 2: LD pruning
echo "Pruning SNPs in LD..."
plink --bfile temp_data \
      --indep-pairwise $WINDOW_SIZE $STEP_SIZE $LD_THRESHOLD \
      --allow-extra-chr \
      --out ld_pruning

# Step 3: Extract pruned SNPs
echo "Extracting pruned SNPs..."
plink --bfile temp_data \
      --extract ld_pruning.prune.in \
      --recode vcf-iid bgz \
      --allow-extra-chr \
      --out ld_pruned.snps

# Step 4: Index the output
echo "Indexing output VCF..."
bcftools index "$LD_FILTERED"

# Check results
TOTAL=$(bcftools view -H "$LD_FILTERED" | wc -l)
echo ""
echo "LD pruning complete!"
echo "Total SNPs: $TOTAL"

# Clean up
rm -f temp_data.* ld_pruning.*


# Calculate and save statistics
STATS_FILE="${LD_FILTERED%.vcf.gz}.stats.txt"
echo "Generating statistics file: $STATS_FILE"
bcftools stats "$LD_FILTERED" > "$STATS_FILE"
echo "Statistics saved to $STATS_FILE"
echo ""

echo "Done!"
}

# Usage: ./subsample_snps.sh <input_vcf> <output_vcf> <number_of_snps>

subsample_snps() {
    local INPUT_VCF=$1
    local OUTPUT_VCF=$2
    local N_SNPS=$3
    
    echo "Subsampling $N_SNPS random SNPs from $INPUT_VCF..."
    
    # Count total SNPs
    TOTAL=$(bcftools view -H "$INPUT_VCF" | wc -l)
    echo "Total SNPs available: $TOTAL"
    
    if [ "$TOTAL" -lt "$N_SNPS" ]; then
        echo "ERROR: Requested $N_SNPS SNPs but only $TOTAL available!"
        return 1
    fi
    
    # Get header
    bcftools view -h "$INPUT_VCF" > temp_header.vcf
    
    # Randomly sample N lines from the variants
    bcftools view -H "$INPUT_VCF" | shuf -n "$N_SNPS" | sort -k1,1 -k2,2n >> temp_header.vcf
    
    # Compress and index
    bgzip -c temp_header.vcf > "$OUTPUT_VCF"
    bcftools index "$OUTPUT_VCF"
    
    # Clean up
    rm temp_header.vcf
    
    # Verify
    FINAL=$(bcftools view -H "$OUTPUT_VCF" | wc -l)
    echo "Created $OUTPUT_VCF with $FINAL SNPs"

    # Calculate and save statistics
    STATS_FILE="${OUTPUT_VCF%.vcf.gz}.stats.txt"
    echo "Generating statistics file: $STATS_FILE"
    bcftools stats "$OUTPUT_VCF" > "$STATS_FILE"
    echo "Statistics saved to $STATS_FILE"
    echo ""
}

# RUN PIPELINE

filter_ld "$IN" "$LD_FILTERED"
subsample_snps "$LD_FILTERED" "$BIG" "50000"
subsample_snps "$LD_FILTERED" "$SMALL" "5000"
