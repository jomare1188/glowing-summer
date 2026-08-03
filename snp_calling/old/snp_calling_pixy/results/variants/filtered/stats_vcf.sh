#!/bin/bash
# Advanced VCF Statistics Extraction with Population Genetics Metrics
# This script uses bcftools and vcftools for comprehensive analysis

VCF_FILE="$1"
OUTPUT_PREFIX="${2:-vcf_stats}"

if [ -z "$VCF_FILE" ]; then
    echo "Usage: $0 <vcf_file> [output_prefix]"
    echo ""
    echo "This script generates comprehensive VCF statistics including:"
    echo "  - Basic variant counts and types"
    echo "  - Quality metrics"
    echo "  - Allele frequency distribution"
    echo "  - Missing data statistics"
    echo "  - Hardy-Weinberg equilibrium tests"
    echo "  - Depth and quality distributions"
    echo ""
    echo "Requirements: bcftools, vcftools (optional for extended stats)"
    exit 1
fi

echo "========================================================================"
echo "VCF COMPREHENSIVE STATISTICS ANALYSIS"
echo "========================================================================"
echo "Input VCF: $VCF_FILE"
echo "Output prefix: $OUTPUT_PREFIX"
echo ""

# Check if bcftools is available
if ! command -v bcftools &> /dev/null; then
    echo "ERROR: bcftools not found. Please install bcftools."
    echo "Install with: conda install -c bioconda bcftools"
    exit 1
fi

echo "1. Extracting basic variant statistics..."
bcftools stats "$VCF_FILE" > "${OUTPUT_PREFIX}.bcftools_stats.txt"

echo "2. Counting variants per chromosome..."
bcftools index -s "$VCF_FILE" 2>/dev/null || echo "Note: VCF not indexed, some operations may be slower"
bcftools query -f '%CHROM\n' "$VCF_FILE" | sort | uniq -c > "${OUTPUT_PREFIX}.variants_per_chromosome.txt"

echo "3. Extracting variant type distribution..."
bcftools query -f '%TYPE\n' "$VCF_FILE" | sort | uniq -c > "${OUTPUT_PREFIX}.variant_types.txt"

echo "4. Calculating allele frequency distribution..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$VCF_FILE" 2>/dev/null > "${OUTPUT_PREFIX}.allele_frequencies.txt" || \
bcftools +fill-tags "$VCF_FILE" -- -t AF,AN,AC | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' > "${OUTPUT_PREFIX}.allele_frequencies.txt"

echo "5. Extracting quality score distribution..."
bcftools query -f '%QUAL\n' "$VCF_FILE" | grep -v '\.' > "${OUTPUT_PREFIX}.quality_scores.txt"

echo "6. Sample-level statistics..."
bcftools query -l "$VCF_FILE" > "${OUTPUT_PREFIX}.sample_list.txt"
NUM_SAMPLES=$(bcftools query -l "$VCF_FILE" | wc -l)
echo "Number of samples: $NUM_SAMPLES" > "${OUTPUT_PREFIX}.sample_count.txt"

# If vcftools is available, run additional analyses
if command -v vcftools &> /dev/null; then
    echo "7. Running vcftools analyses..."
    
    echo "   - Missing data per individual..."
    vcftools --gzvcf "$VCF_FILE" --missing-indv --out "${OUTPUT_PREFIX}" 2>&1 | grep -v "^After filtering"
    
    echo "   - Missing data per site..."
    vcftools --gzvcf "$VCF_FILE" --missing-site --out "${OUTPUT_PREFIX}" 2>&1 | grep -v "^After filtering"
    
    echo "   - Site depth statistics..."
    vcftools --gzvcf "$VCF_FILE" --site-mean-depth --out "${OUTPUT_PREFIX}" 2>&1 | grep -v "^After filtering"
    
    echo "   - Allele frequency..."
    vcftools --gzvcf "$VCF_FILE" --freq2 --out "${OUTPUT_PREFIX}" 2>&1 | grep -v "^After filtering"
    
    echo "   - Site quality scores..."
    vcftools --gzvcf "$VCF_FILE" --site-quality --out "${OUTPUT_PREFIX}" 2>&1 | grep -v "^After filtering"
    
    # Hardy-Weinberg if applicable
    echo "   - Hardy-Weinberg equilibrium..."
    vcftools --gzvcf "$VCF_FILE" --hardy --out "${OUTPUT_PREFIX}" 2>&1 | grep -v "^After filtering"
    
    echo "   - Heterozygosity per individual..."
    vcftools --gzvcf "$VCF_FILE" --het --out "${OUTPUT_PREFIX}" 2>&1 | grep -v "^After filtering"
else
    echo "7. vcftools not found - skipping extended statistics"
    echo "   Install with: conda install -c bioconda vcftools"
fi

echo ""
echo "========================================================================"
echo "GENERATING SUMMARY REPORT"
echo "========================================================================"

# Parse bcftools stats for key metrics
echo "Extracting key statistics from bcftools stats..."

cat << 'EOF' > "${OUTPUT_PREFIX}.summary_report.txt"
======================================================================
VCF STATISTICS SUMMARY REPORT
======================================================================

EOF

echo "File: $VCF_FILE" >> "${OUTPUT_PREFIX}.summary_report.txt"
echo "Analysis date: $(date)" >> "${OUTPUT_PREFIX}.summary_report.txt"
echo "" >> "${OUTPUT_PREFIX}.summary_report.txt"

# Extract number of samples
echo "Number of samples: $NUM_SAMPLES" >> "${OUTPUT_PREFIX}.summary_report.txt"
echo "" >> "${OUTPUT_PREFIX}.summary_report.txt"

# Extract SNP summary from bcftools stats
echo "VARIANT SUMMARY:" >> "${OUTPUT_PREFIX}.summary_report.txt"
echo "----------------------------------------------------------------------" >> "${OUTPUT_PREFIX}.summary_report.txt"
grep "^SN" "${OUTPUT_PREFIX}.bcftools_stats.txt" | cut -f 2- >> "${OUTPUT_PREFIX}.summary_report.txt"
echo "" >> "${OUTPUT_PREFIX}.summary_report.txt"

# Extract TSTV (transition/transversion) ratio
echo "TRANSITION/TRANSVERSION RATIO:" >> "${OUTPUT_PREFIX}.summary_report.txt"
echo "----------------------------------------------------------------------" >> "${OUTPUT_PREFIX}.summary_report_report.txt"
grep "^TSTV" "${OUTPUT_PREFIX}.bcftools_stats.txt" | cut -f 2- >> "${OUTPUT_PREFIX}.summary_report.txt"
echo "" >> "${OUTPUT_PREFIX}.summary_report.txt"

# Variant types distribution
echo "VARIANT TYPES DISTRIBUTION:" >> "${OUTPUT_PREFIX}.summary_report.txt"
echo "----------------------------------------------------------------------" >> "${OUTPUT_PREFIX}.summary_report.txt"
cat "${OUTPUT_PREFIX}.variant_types.txt" >> "${OUTPUT_PREFIX}.summary_report.txt"
echo "" >> "${OUTPUT_PREFIX}.summary_report.txt"

# Variants per chromosome
echo "VARIANTS PER CHROMOSOME:" >> "${OUTPUT_PREFIX}.summary_report.txt"
echo "----------------------------------------------------------------------" >> "${OUTPUT_PREFIX}.summary_report.txt"
cat "${OUTPUT_PREFIX}.variants_per_chromosome.txt" >> "${OUTPUT_PREFIX}.summary_report.txt"
echo "" >> "${OUTPUT_PREFIX}.summary_report.txt"

echo "======================================================================" >> "${OUTPUT_PREFIX}.summary_report.txt"
echo "FILES GENERATED:" >> "${OUTPUT_PREFIX}.summary_report.txt"
echo "======================================================================" >> "${OUTPUT_PREFIX}.summary_report.txt"
ls -lh "${OUTPUT_PREFIX}".* >> "${OUTPUT_PREFIX}.summary_report.txt"

echo ""
echo "Analysis complete!"
echo ""
echo "Key output files:"
echo "  - ${OUTPUT_PREFIX}.summary_report.txt (main summary)"
echo "  - ${OUTPUT_PREFIX}.bcftools_stats.txt (detailed bcftools statistics)"
echo "  - ${OUTPUT_PREFIX}.allele_frequencies.txt (AF for each variant)"
echo "  - ${OUTPUT_PREFIX}.variants_per_chromosome.txt (counts by chromosome)"
if command -v vcftools &> /dev/null; then
    echo "  - ${OUTPUT_PREFIX}.imiss (missing data per individual)"
    echo "  - ${OUTPUT_PREFIX}.lmiss (missing data per site)"
    echo "  - ${OUTPUT_PREFIX}.het (heterozygosity per individual)"
fi
echo ""
echo "View summary with: cat ${OUTPUT_PREFIX}.summary_report.txt"
echo ""
