#!/usr/bin/env bash
# =============================================================================
# ROH ANALYSIS PIPELINE — Callicarpa sp.
# Individual-level inbreeding audit for conservation genomics
# Genome size: 506,362,408 bp
# =============================================================================
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION — edit these before running
# ─────────────────────────────────────────────────────────────────────────────
INPUT_VCF="../cohort.snps.pass.vcf.gz"     # Your pre-subsampling filtered VCF
GENOME_SIZE=506362408                    # Callicarpa genome size in bp
OUTDIR="roh_analysis"
LOG="${OUTDIR}/pipeline.log"

# ─────────────────────────────────────────────────────────────────────────────
# SETUP
# ─────────────────────────────────────────────────────────────────────────────
mkdir -p "${OUTDIR}"/{density_check,maf_comparison,roh_short,roh_medium,roh_long,froh_results}

exec > >(tee -a "${LOG}") 2>&1
echo "============================================================"
echo " ROH Pipeline started: $(date)"
echo " Input VCF: ${INPUT_VCF}"
echo " Genome size: ${GENOME_SIZE} bp"
echo "============================================================"


# =============================================================================
# STEP 1: PREPARE FILTERED VCF (two MAF thresholds, NO HWE filter)
# =============================================================================
echo ""
echo ">>> STEP 1: VCF filtering (MAF 0.05 and MAF 0.10)"
echo "    NOTE: HWE filter intentionally excluded — it removes inbreeding signal"

for MAF in 0.05 0.10; do
    TAG="maf${MAF/./}"   # e.g. maf005, maf010
    OUT_VCF="${OUTDIR}/maf_comparison/for_roh_dense_${TAG}.vcf.gz"

    BIALLELIC_TMP="${OUTDIR}/maf_comparison/biallelic_tmp_${TAG}.vcf.gz"

    # --- Sub-step A: extract biallelic SNPs only into a real intermediate file ---
    # Piping bcftools->vcftools is unreliable: vcftools silently falls back to
    # reading /dev/stdin which can fail, leaving multiallelic sites in the output.
    echo "  -> [A] Extracting biallelic SNPs only (bcftools) ..."
    bcftools view \
        --min-alleles 2 \
        --max-alleles 2 \
        --type snps \
        -O z -o "${BIALLELIC_TMP}" \
        "${INPUT_VCF}"
    tabix "${BIALLELIC_TMP}"

    N_BIALLELIC=$(bcftools stats "${BIALLELIC_TMP}" | grep "^SN.*number of SNPs" | awk '{print $NF}')
    echo "     Biallelic SNPs before MAF/missingness filter: ${N_BIALLELIC}"

    # --- Sub-step B: apply call-rate + MAF filter with vcftools ---
    echo "  -> [B] Applying MAF=${MAF} and call-rate filter (vcftools) ..."
    vcftools \
        --gzvcf "${BIALLELIC_TMP}" \
        --max-missing 0.9 \
        --maf "${MAF}" \
        --recode --stdout \
    | bgzip > "${OUT_VCF}"

    tabix "${OUT_VCF}"
    rm -f "${BIALLELIC_TMP}" "${BIALLELIC_TMP}.tbi"   # clean up temp

    N_SNPS=$(bcftools stats "${OUT_VCF}" | grep "^SN.*number of SNPs" | awk '{print $NF}')
    echo "     MAF=${MAF}: ${N_SNPS} SNPs retained in final filtered VCF"
done


# =============================================================================
# STEP 2: SNP DENSITY CHECK
# =============================================================================
echo ""
echo ">>> STEP 2: SNP density diagnostics"

for MAF in 0.05 0.10; do
    TAG="maf${MAF/./}"
    VCF="${OUTDIR}/maf_comparison/for_roh_dense_${TAG}.vcf.gz"
    DENSITY_OUT="${OUTDIR}/density_check/density_${TAG}.txt"

    echo "  -> Density check for MAF=${MAF} ..."

    # SNPs per chromosome + average spacing
    bcftools query -f '%CHROM\t%POS\n' "${VCF}" \
    | awk '
        {
            chrom = $1; pos = $2
            count[chrom]++
            if (!(chrom in first)) first[chrom] = pos
            last[chrom] = pos
        }
        END {
            printf "%-20s %10s %15s %15s %15s\n",
                "CHROM","N_SNPS","FIRST_POS","LAST_POS","AVG_SPACING_KB"
            for (c in count) {
                span = last[c] - first[c]
                spacing = (count[c] > 1) ? (span / (count[c]-1)) / 1000 : 0
                printf "%-20s %10d %15d %15d %15.2f\n",
                    c, count[c], first[c], last[c], spacing
            }
        }
    ' | sort -k1,1 > "${DENSITY_OUT}"

    echo "     Results saved to ${DENSITY_OUT}"
    cat "${DENSITY_OUT}"
    echo ""
done


# =============================================================================
# STEP 3: CONVERT TO PLINK FORMAT (both MAF VCFs)
# =============================================================================
echo ""
echo ">>> STEP 3: Converting to PLINK binary format"

for MAF in 0.05 0.10; do
    TAG="maf${MAF/./}"
    VCF="${OUTDIR}/maf_comparison/for_roh_dense_${TAG}.vcf.gz"
    PLINK_PREFIX="${OUTDIR}/maf_comparison/callicarpa_${TAG}"

    echo "  -> plink2 conversion for MAF=${MAF} ..."
    plink2 \
        --vcf "${VCF}" \
        --make-bed \
        --out "${PLINK_PREFIX}" \
        --allow-extra-chr \
        --max-alleles 2 \
        --set-missing-var-ids '@:#' \
        --new-id-max-allele-len 50

    echo "     PLINK files: ${PLINK_PREFIX}.{bed,bim,fam}"
done


# =============================================================================
# STEP 4: ROH DETECTION — THREE LENGTH CLASSES × TWO MAF THRESHOLDS
# =============================================================================
# Length classes explained:
#   SHORT  (>100 kb)  → ancient/background inbreeding, many generations ago
#   MEDIUM (>500 kb)  → moderate; several generations back (your original target)
#   LONG   (>1500 kb) → recent/close inbreeding, 1–3 generations

echo ""
echo ">>> STEP 4: ROH detection — 3 length classes × 2 MAF thresholds"

declare -A ROH_KB=(     [short]=100    [medium]=500    [long]=1500  )
declare -A ROH_SNP=(    [short]=25     [medium]=50     [long]=100   )
declare -A ROH_GAP=(    [short]=1000   [medium]=1000   [long]=500   )
declare -A ROH_WINSNP=( [short]=25     [medium]=50     [long]=100   )

for CLASS in short medium long; do
    KB=${ROH_KB[$CLASS]}
    SNP=${ROH_SNP[$CLASS]}
    GAP=${ROH_GAP[$CLASS]}
    WINSNP=${ROH_WINSNP[$CLASS]}

    for MAF in 0.05 0.10; do
        TAG="maf${MAF/./}"
        PLINK_IN="${OUTDIR}/maf_comparison/callicarpa_${TAG}"
        OUT_PREFIX="${OUTDIR}/roh_${CLASS}/callicarpa_${CLASS}_${TAG}"

        echo "  -> ROH ${CLASS} (>=${KB}kb, snp>=${SNP}, gap<=${GAP}kb) | MAF=${MAF}"

        plink \
            --bfile "${PLINK_IN}" \
            --homozyg \
            --homozyg-kb      "${KB}" \
            --homozyg-snp     "${SNP}" \
            --homozyg-density 50 \
            --homozyg-gap     "${GAP}" \
            --homozyg-window-snp     "${WINSNP}" \
            --homozyg-window-het     1 \
            --homozyg-window-missing 5 \
            --allow-extra-chr \
            --out "${OUT_PREFIX}" \
        && echo "     Done -> ${OUT_PREFIX}.hom.indiv" \
        || echo "     WARNING: PLINK exited with error for ${CLASS}/${TAG}"

    done
done


# =============================================================================
# STEP 5: F_ROH CALCULATION IN R
# =============================================================================
echo ""
echo ">>> STEP 5: Computing F_ROH per individual"

R --no-save --quiet << RSCRIPT
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

GENOME_SIZE <- ${GENOME_SIZE}
OUTDIR      <- "${OUTDIR}"

classes <- c("short", "medium", "long")
mafs    <- c("maf005", "maf010")
results <- list()

for (cls in classes) {
  for (maf in mafs) {
    f <- file.path(OUTDIR, paste0("roh_", cls),
                   paste0("callicarpa_", cls, "_", maf, ".hom.indiv"))
    if (!file.exists(f)) {
      message("Skipping (not found): ", f); next
    }
    df <- read_table(f, show_col_types = FALSE)
    df\$CLASS   <- cls
    df\$MAF_TAG <- maf
    df\$FROH    <- (df\$KB * 1000) / GENOME_SIZE
    results[[paste(cls, maf)]] <- df
  }
}

all <- bind_rows(results)

# ── Save master table ──────────────────────────────────────────────────────
write_csv(all, file.path(OUTDIR, "froh_results", "froh_all_classes.csv"))
message("Saved: froh_all_classes.csv")

# ── Per-individual summary across length classes (MAF 0.05, as recommended) ─
summary_df <- all %>%
  filter(MAF_TAG == "maf005") %>%
  select(IID, CLASS, FROH, NSEG, KB) %>%
  pivot_wider(
    names_from  = CLASS,
    values_from = c(FROH, NSEG, KB),
    names_sep   = "_"
  ) %>%
  arrange(desc(FROH_medium))

write_csv(summary_df, file.path(OUTDIR, "froh_results", "froh_per_individual_summary.csv"))
message("Saved: froh_per_individual_summary.csv")
print(summary_df)

# ── MAF sensitivity check (medium ROH class) ─────────────────────────────
maf_compare <- all %>%
  filter(CLASS == "medium") %>%
  select(IID, MAF_TAG, FROH) %>%
  pivot_wider(names_from = MAF_TAG, values_from = FROH, names_prefix = "FROH_") %>%
  mutate(delta_FROH = abs(FROH_maf005 - FROH_maf010))

write_csv(maf_compare, file.path(OUTDIR, "froh_results", "froh_maf_sensitivity.csv"))
message("MAF sensitivity check (medium ROH):")
print(maf_compare)

# ── PLOTS ──────────────────────────────────────────────────────────────────

# 1. F_ROH by individual and ROH class (MAF=0.05)
p1 <- all %>%
  filter(MAF_TAG == "maf005") %>%
  mutate(CLASS = factor(CLASS, levels = c("short","medium","long"))) %>%
  ggplot(aes(x = reorder(IID, -FROH), y = FROH, fill = CLASS)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(short="#a8c7fa", medium="#4a90d9", long="#1a3a6b"),
                    labels = c("Short (>100kb)","Medium (>500kb)","Long (>1500kb)")) +
  labs(
    title    = "F_ROH per Individual — Callicarpa sp.",
    subtitle = paste0("Genome size: ", format(GENOME_SIZE, big.mark=","), " bp | MAF ≥ 0.05"),
    x        = "Individual",
    y        = expression(F[ROH]),
    fill     = "ROH Class"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTDIR, "froh_results", "froh_by_individual.pdf"),
       p1, width = 12, height = 6)
ggsave(file.path(OUTDIR, "froh_results", "froh_by_individual.png"),
       p1, width = 12, height = 6, dpi = 150)
message("Saved plot: froh_by_individual.pdf / .png")

# 2. F_ROH distribution across individuals (boxplot by class)
p2 <- all %>%
  filter(MAF_TAG == "maf005") %>%
  mutate(CLASS = factor(CLASS, levels = c("short","medium","long"))) %>%
  ggplot(aes(x = CLASS, y = FROH, fill = CLASS)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c(short="#a8c7fa", medium="#4a90d9", long="#1a3a6b")) +
  scale_x_discrete(labels = c("Short\n(>100kb)","Medium\n(>500kb)","Long\n(>1500kb)")) +
  labs(
    title = "F_ROH Distribution by ROH Length Class",
    x = "ROH Length Class", y = expression(F[ROH])
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(OUTDIR, "froh_results", "froh_distribution_by_class.pdf"),
       p2, width = 7, height = 5)
ggsave(file.path(OUTDIR, "froh_results", "froh_distribution_by_class.png"),
       p2, width = 7, height = 5, dpi = 150)
message("Saved plot: froh_distribution_by_class.pdf / .png")

# 3. MAF sensitivity scatter (medium ROH)
p3 <- maf_compare %>%
  ggplot(aes(x = FROH_maf010, y = FROH_maf005, label = IID)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, color = "#4a90d9") +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
  labs(
    title    = "MAF Sensitivity Check — Medium ROH (>500kb)",
    subtitle = "Points above the line: MAF 0.05 detects more inbreeding than MAF 0.10",
    x        = expression(F[ROH] ~ "(MAF ≥ 0.10)"),
    y        = expression(F[ROH] ~ "(MAF ≥ 0.05)")
  ) +
  theme_bw(base_size = 12)

# ggrepel may not be installed — fallback
tryCatch({
  ggsave(file.path(OUTDIR, "froh_results", "froh_maf_sensitivity_scatter.pdf"),
         p3, width = 7, height = 6)
  ggsave(file.path(OUTDIR, "froh_results", "froh_maf_sensitivity_scatter.png"),
         p3, width = 7, height = 6, dpi = 150)
}, error = function(e) {
    p3b <- p3 + geom_text(size = 2.5, vjust = -0.5)
    ggsave(file.path(OUTDIR, "froh_results", "froh_maf_sensitivity_scatter.pdf"),
           p3b, width = 7, height = 6)
    ggsave(file.path(OUTDIR, "froh_results", "froh_maf_sensitivity_scatter.png"),
           p3b, width = 7, height = 6, dpi = 150)
})
message("Saved plot: froh_maf_sensitivity_scatter.pdf / .png")

message("")
message("============================================================")
message(" F_ROH SUMMARY (MAF=0.05, medium ROH class >500kb)")
message("============================================================")
medium_summary <- all %>%
  filter(MAF_TAG == "maf005", CLASS == "medium") %>%
  summarise(
    N_individuals  = n(),
    mean_FROH      = round(mean(FROH), 4),
    median_FROH    = round(median(FROH), 4),
    max_FROH       = round(max(FROH), 4),
    min_FROH       = round(min(FROH), 4),
    pct_FROH_gt0.1 = round(mean(FROH > 0.1) * 100, 1)
  )
print(medium_summary)
message("  Note: F_ROH > 0.125 is roughly equivalent to first-cousin inbreeding")
message("  Note: F_ROH > 0.25  indicates parent-offspring or full-sib mating")

RSCRIPT

echo ""
echo "============================================================"
echo " Pipeline complete: $(date)"
echo " All outputs in: ${OUTDIR}/"
echo ""
echo " Key files:"
echo "   ${OUTDIR}/density_check/          — SNP density per chromosome"
echo "   ${OUTDIR}/roh_*/                  — PLINK .hom and .hom.indiv files"
echo "   ${OUTDIR}/froh_results/froh_all_classes.csv"
echo "   ${OUTDIR}/froh_results/froh_per_individual_summary.csv"
echo "   ${OUTDIR}/froh_results/froh_maf_sensitivity.csv"
echo "   ${OUTDIR}/froh_results/froh_by_individual.pdf"
echo "   ${OUTDIR}/froh_results/froh_distribution_by_class.pdf"
echo "   ${OUTDIR}/froh_results/froh_maf_sensitivity_scatter.pdf"
echo "============================================================"
