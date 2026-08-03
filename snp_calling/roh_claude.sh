#!/usr/bin/env bash
################################################################################
# ROH ANALYSIS — Callicarpa americana
# Individual-level inbreeding audit for conservation genomics
#
# Consumes the depth-equalised datasets built by snp_datasets.sh. The VCF
# filtering and PLINK conversion that this script used to do itself are now
# done upstream (snp_datasets.sh Phase 4), so this script starts from the
# already-prepared PLINK binaries and does only what is specific to ROH.
#
# Inputs, per dataset:
#   results/<tag>/roh/callicarpa_maf{005,010}.{bed,bim,fam}
#   results/<tag>/roh/for_roh_dense_maf{005,010}.vcf.gz   (density check only)
#
# Usage:
#   bash roh_claude.sh              # both datasets, then the cross comparison
#   bash roh_claude.sh d09          # one dataset only
#   FORCE=1 bash roh_claude.sh      # recompute even where outputs exist
#
# Note on comparing datasets: a sample present in both was called twice, from
# two independently downsampled BAMs. d09 and d20 are therefore NOT nested
# subsets of one call set, and F_ROH may legitimately differ between them.
# That difference is the point of the cross-dataset comparison in step 4.
################################################################################
set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────
PROJECT_DIR="/dados04/jorge/CALLICARPA_2"
WORK_DIR="${PROJECT_DIR}/snp_calling"
OUT_DIR="${WORK_DIR}/results"
GENOME_SIZE=506362408

DATASETS=("d09:9" "d20:20")
ROH_MAFS=(0.05 0.10)
PRIMARY_MAF_TAG="maf005"        # reported as the primary result
FORCE=${FORCE:-0}

# ROH length classes. The long class uses a stricter gap (500 kb vs 1000 kb)
# deliberately: merging two genuinely distinct recent segments would inflate
# the very signal the long class exists to measure.
declare -A ROH_KB=(     [short]=100    [medium]=500    [long]=1500  )
declare -A ROH_SNP=(    [short]=25     [medium]=50     [long]=100   )
declare -A ROH_GAP=(    [short]=1000   [medium]=1000   [long]=500   )
declare -A ROH_WINSNP=( [short]=25     [medium]=50     [long]=100   )
ROH_CLASSES=(short medium long)

# plink here must be 1.9: --homozyg does not exist in PLINK 2. popgen_tools is
# prepended so `plink` resolves to 1.9 even though SNP_call ships plink2.
CONDA_BASE="/home/genomics/miniconda3"
export PATH="${CONDA_BASE}/envs/popgen_tools/bin:${CONDA_BASE}/envs/SNP_call/bin:${PATH}"

# SNP_call's R is 3.5.0 and has only dplyr; the plotting/tidy stack lives in
# R_popstat_jorge (R 4.4.3), which is the only env carrying all of
# dplyr/tidyr/readr/ggplot2/ggrepel.
RSCRIPT_BIN="${RSCRIPT_BIN:-${CONDA_BASE}/envs/R_popstat_jorge/bin/Rscript}"

LOG_DIR="${WORK_DIR}/logs"
mkdir -p "${LOG_DIR}"
MAIN_LOG="${LOG_DIR}/roh_analysis.log"

# ─────────────────────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

_emit() { echo -e "$1[$(date '+%F %T')] $2${NC}" | tee -a "${MAIN_LOG}"; }
log()   { _emit "${BLUE}"   "$*"; }
info()  { _emit ""          "  $*"; }
ok()    { _emit "${GREEN}"  "  $*"; }
warn()  { _emit "${YELLOW}" "  WARNING: $*"; }
die()   { _emit "${RED}"    "ERROR: $*"; exit 1; }

maf_tag() { echo "maf${1/./}"; }

# ─────────────────────────────────────────────────────────────────────────────
# PREFLIGHT
# ─────────────────────────────────────────────────────────────────────────────
preflight() {
    log "Preflight"

    for t in plink bcftools; do
        command -v "$t" >/dev/null || die "$t not found in PATH"
    done

    # PLINK 2 silently lacks --homozyg, so verify we really have 1.9.
    local pv
    pv=$(plink --version 2>&1 | head -1)
    [[ "$pv" == *"v1.9"* ]] || die "plink is not v1.9 (--homozyg unavailable): ${pv}"
    info "plink: ${pv}"

    [[ -x "${RSCRIPT_BIN}" ]] || die "Rscript not found: ${RSCRIPT_BIN}"
    local missing
    missing=$("${RSCRIPT_BIN}" -e 'cat(paste(setdiff(c("dplyr","tidyr","readr","ggplot2"), rownames(installed.packages())), collapse=" "))' 2>/dev/null)
    [[ -z "${missing}" ]] || die "R packages missing from ${RSCRIPT_BIN}: ${missing}"
    info "Rscript: ${RSCRIPT_BIN}"

    for spec in "${TARGETS[@]}"; do
        local tag="${spec%%:*}"
        for maf in "${ROH_MAFS[@]}"; do
            local mt; mt=$(maf_tag "$maf")
            local pfx="${OUT_DIR}/${tag}/roh/callicarpa_${mt}"
            for ext in bed bim fam; do
                [[ -s "${pfx}.${ext}" ]] || die "missing ${pfx}.${ext} — run: bash snp_datasets.sh downstream"
            done
        done
        info "${tag}: PLINK inputs present ($(wc -l < "${OUT_DIR}/${tag}/roh/callicarpa_${PRIMARY_MAF_TAG}.fam") samples)"
    done
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: SNP DENSITY DIAGNOSTICS
# ─────────────────────────────────────────────────────────────────────────────
# ROH parameters are only defensible if marker density supports them: a
# --homozyg-kb 100 window needs enough SNPs in 100 kb to distinguish a real
# tract from a run of missing data.
density_check() {
    local tag=$1
    local rdir="${OUT_DIR}/${tag}/roh_analysis"
    mkdir -p "${rdir}/density_check"

    for maf in "${ROH_MAFS[@]}"; do
        local mt; mt=$(maf_tag "$maf")
        local vcf="${OUT_DIR}/${tag}/roh/for_roh_dense_${mt}.vcf.gz"
        local out="${rdir}/density_check/density_${mt}.txt"

        if [[ -s "${out}" && "${FORCE}" != "1" ]]; then
            info "${tag}/${mt}: density check exists, skipping"
            continue
        fi
        [[ -s "${vcf}" ]] || { warn "${tag}: ${vcf} absent, skipping density check"; continue; }

        bcftools query -f '%CHROM\t%POS\n' "${vcf}" \
        | awk '
            { c=$1; p=$2; n[c]++; if (!(c in first)) first[c]=p; last[c]=p }
            END {
                printf "%-24s %10s %14s %14s %16s\n", "CHROM","N_SNPS","FIRST_POS","LAST_POS","AVG_SPACING_KB"
                for (c in n) {
                    sp = (n[c] > 1) ? ((last[c]-first[c]) / (n[c]-1)) / 1000 : 0
                    printf "%-24s %10d %14d %14d %16.2f\n", c, n[c], first[c], last[c], sp
                }
            }' \
        | { read -r hdr; echo "$hdr"; sort -k2,2nr; } > "${out}"

        local nscaf med
        nscaf=$(( $(wc -l < "${out}") - 1 ))
        med=$(awk 'NR>1 {print $5}' "${out}" | sort -n | awk '{a[NR]=$1} END{print (NR%2 ? a[(NR+1)/2] : (a[NR/2]+a[NR/2+1])/2)}')
        ok "${tag}/${mt}: ${nscaf} scaffolds, median spacing ${med} kb"
    done
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: ROH DETECTION — 3 LENGTH CLASSES x 2 MAF THRESHOLDS
# ─────────────────────────────────────────────────────────────────────────────
detect_roh() {
    local tag=$1
    local rdir="${OUT_DIR}/${tag}/roh_analysis"

    for cls in "${ROH_CLASSES[@]}"; do
        mkdir -p "${rdir}/roh_${cls}"
        local kb=${ROH_KB[$cls]} snp=${ROH_SNP[$cls]} gap=${ROH_GAP[$cls]} win=${ROH_WINSNP[$cls]}

        for maf in "${ROH_MAFS[@]}"; do
            local mt; mt=$(maf_tag "$maf")
            local pin="${OUT_DIR}/${tag}/roh/callicarpa_${mt}"
            local pout="${rdir}/roh_${cls}/callicarpa_${cls}_${mt}"

            if [[ -s "${pout}.hom.indiv" && "${FORCE}" != "1" ]]; then
                info "${tag}/${cls}/${mt}: exists, skipping"
                continue
            fi

            info "${tag}/${cls}/${mt}: >=${kb}kb, >=${snp} SNPs, gap <=${gap}kb"
            if plink --bfile "${pin}" --homozyg \
                     --homozyg-kb "${kb}" --homozyg-snp "${snp}" \
                     --homozyg-density 50 --homozyg-gap "${gap}" \
                     --homozyg-window-snp "${win}" \
                     --homozyg-window-het 1 --homozyg-window-missing 5 \
                     --allow-extra-chr --out "${pout}" \
                     &> "${LOG_DIR}/${tag}_roh_${cls}_${mt}.log"
            then
                local nseg
                nseg=$(( $(wc -l < "${pout}.hom") - 1 ))
                ok "${tag}/${cls}/${mt}: ${nseg} segments"
            else
                warn "${tag}/${cls}/${mt}: plink failed, see ${LOG_DIR}/${tag}_roh_${cls}_${mt}.log"
            fi
        done
    done
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 & 4: F_ROH AND CROSS-DATASET COMPARISON (R)
# ─────────────────────────────────────────────────────────────────────────────
# The R source is written out verbatim (quoted heredoc) and driven by argv, so
# no shell expansion happens inside it. The previous version interpolated the
# heredoc and needed \$ escapes on every dplyr column reference.
write_r_script() {
    cat > "$1" <<'RSCRIPT'
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

args        <- commandArgs(trailingOnly = TRUE)
OUT_DIR     <- args[1]
GENOME_SIZE <- as.numeric(args[2])
PRIMARY_MAF <- args[3]
tags        <- args[-(1:3)]

classes <- c("short", "medium", "long")
mafs    <- c("maf005", "maf010")
pal     <- c(short = "#a8c7fa", medium = "#4a90d9", long = "#1a3a6b")
lab_cls <- c("Short (>100kb)", "Medium (>500kb)", "Long (>1500kb)")

has_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (!has_repel) message("Note: ggrepel not installed; using plain labels")

# PDF for print, PNG so the figures can be embedded in the README — GitHub
# renders neither PDF nor SVG inline in markdown.
save_plot <- function(stem, p, w, h) {
  ggsave(paste0(stem, ".pdf"), p, width = w, height = h)
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 150, bg = "white")
}

read_dataset <- function(tag) {
  rows <- list()
  for (cls in classes) for (mf in mafs) {
    f <- file.path(OUT_DIR, tag, "roh_analysis", paste0("roh_", cls),
                   paste0("callicarpa_", cls, "_", mf, ".hom.indiv"))
    if (!file.exists(f)) { message("  skipping (absent): ", f); next }
    # base read.table, not readr: PLINK pads .hom.indiv rows with a trailing
    # space, which read_table tokenises as a 7th column against a 6-name
    # header and reports as a parsing failure on every row.
    d <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
    d$CLASS <- cls; d$MAF_TAG <- mf; d$DATASET <- tag
    d$FROH  <- (d$KB * 1000) / GENOME_SIZE
    rows[[paste(cls, mf)]] <- d
  }
  if (length(rows) == 0) return(NULL)
  bind_rows(rows)
}

all_list <- list()

for (tag in tags) {
  message("\n=== F_ROH: ", tag, " ===")
  a <- read_dataset(tag)
  if (is.null(a)) { message("  no .hom.indiv files, skipping"); next }
  all_list[[tag]] <- a

  fdir <- file.path(OUT_DIR, tag, "roh_analysis", "froh_results")
  dir.create(fdir, recursive = TRUE, showWarnings = FALSE)

  write_csv(a, file.path(fdir, "froh_all_classes.csv"))

  summary_df <- a %>%
    filter(MAF_TAG == PRIMARY_MAF) %>%
    select(IID, CLASS, FROH, NSEG, KB) %>%
    pivot_wider(names_from = CLASS, values_from = c(FROH, NSEG, KB), names_sep = "_") %>%
    arrange(desc(FROH_medium))
  write_csv(summary_df, file.path(fdir, "froh_per_individual_summary.csv"))
  print(as.data.frame(summary_df))

  maf_cmp <- a %>%
    filter(CLASS == "medium") %>%
    select(IID, MAF_TAG, FROH) %>%
    pivot_wider(names_from = MAF_TAG, values_from = FROH, names_prefix = "FROH_") %>%
    mutate(delta_FROH = abs(FROH_maf005 - FROH_maf010))
  write_csv(maf_cmp, file.path(fdir, "froh_maf_sensitivity.csv"))

  ap <- a %>% filter(MAF_TAG == PRIMARY_MAF) %>%
        mutate(CLASS = factor(CLASS, levels = classes))

  p1 <- ggplot(ap, aes(reorder(IID, -FROH), FROH, fill = CLASS)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = pal, labels = lab_cls) +
    labs(title = paste0("F_ROH per Individual - ", tag),
         subtitle = paste0(nrow(ap) / length(classes), " individuals | MAF >= 0.05 | genome ",
                           format(GENOME_SIZE, big.mark = ","), " bp"),
         x = "Individual", y = expression(F[ROH]), fill = "ROH Class") +
    theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(file.path(fdir, "froh_by_individual"), p1, 12, 6)

  p2 <- ggplot(ap, aes(CLASS, FROH, fill = CLASS)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    geom_jitter(width = 0.1, size = 2, alpha = 0.8) +
    scale_fill_manual(values = pal) +
    scale_x_discrete(labels = c("Short\n(>100kb)", "Medium\n(>500kb)", "Long\n(>1500kb)")) +
    labs(title = paste0("F_ROH Distribution by ROH Length Class - ", tag),
         x = "ROH Length Class", y = expression(F[ROH])) +
    theme_bw(base_size = 12) + theme(legend.position = "none")
  save_plot(file.path(fdir, "froh_distribution_by_class"), p2, 7, 5)

  p3 <- ggplot(maf_cmp, aes(FROH_maf010, FROH_maf005, label = IID)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 3, colour = "#4a90d9") +
    labs(title = paste0("MAF Sensitivity - Medium ROH (>500kb) - ", tag),
         subtitle = "Above the line: MAF 0.05 detects more inbreeding than MAF 0.10",
         x = expression(F[ROH] ~ "(MAF >= 0.10)"), y = expression(F[ROH] ~ "(MAF >= 0.05)")) +
    theme_bw(base_size = 12)
  p3 <- if (has_repel) p3 + ggrepel::geom_text_repel(size = 3, max.overlaps = 20)
        else           p3 + geom_text(size = 2.5, vjust = -0.5)
  save_plot(file.path(fdir, "froh_maf_sensitivity_scatter"), p3, 7, 6)

  med <- a %>% filter(MAF_TAG == PRIMARY_MAF, CLASS == "medium") %>%
    summarise(N = n(), mean_FROH = round(mean(FROH), 4), median_FROH = round(median(FROH), 4),
              min_FROH = round(min(FROH), 4), max_FROH = round(max(FROH), 4),
              pct_gt_0.125 = round(mean(FROH > 0.125) * 100, 1))
  message("Medium-class summary (", PRIMARY_MAF, "):")
  print(as.data.frame(med))
}

# ── Cross-dataset comparison ────────────────────────────────────────────────
# Samples in both datasets were called twice from independently downsampled
# BAMs. Agreement here is evidence that F_ROH is driven by genuine IBD tracts
# and not by depth; systematic disagreement means depth still matters and the
# deeper dataset should be preferred.
if (length(all_list) >= 2) {
  message("\n=== Cross-dataset comparison ===")
  cdir <- file.path(OUT_DIR, "roh_comparison")
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)

  comb <- bind_rows(all_list) %>% filter(MAF_TAG == PRIMARY_MAF)
  write_csv(comb, file.path(cdir, "froh_all_datasets.csv"))

  wide <- comb %>%
    select(IID, DATASET, CLASS, FROH) %>%
    pivot_wider(names_from = DATASET, values_from = FROH, names_prefix = "FROH_")

  dcols <- setdiff(names(wide), c("IID", "CLASS"))
  dcols_tag <- sub("^FROH_", "", dcols)   # "d09", "d20" — for labels and colours
  if (length(dcols) == 2) {
    wide <- wide %>%
      mutate(delta = .data[[dcols[2]]] - .data[[dcols[1]]]) %>%
      filter(!is.na(.data[[dcols[1]]]), !is.na(.data[[dcols[2]]]))
    write_csv(wide, file.path(cdir, "froh_dataset_concordance.csv"))

    message("Shared individuals: ", length(unique(wide$IID)))
    # cor() warns and returns NA when a class is all-zero in both datasets
    # (e.g. no long ROH detected anywhere), so guard on variance first.
    safe_cor <- function(x, y) {
      if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) return(NA_real_)
      round(cor(x, y), 3)
    }
    conc <- wide %>% group_by(CLASS) %>%
      summarise(n = n(),
                r = safe_cor(.data[[dcols[1]]], .data[[dcols[2]]]),
                mean_abs_delta = round(mean(abs(delta)), 4), .groups = "drop")
    print(as.data.frame(conc))

    # ── PLOT 1: the results themselves, both datasets side by side ──────────
    # Ordered by the d20 value where it exists, else d09, so the panel reads
    # as a ranking. Free y per class: short is ~30x medium and long is zero,
    # so a shared scale would flatten everything but the short panel.
    ord <- comb %>% filter(CLASS == "short") %>%
      group_by(IID) %>%
      summarise(key = ifelse(any(DATASET == dcols_tag[2]),
                             FROH[DATASET == dcols_tag[2]][1], FROH[1] / 2), .groups = "drop") %>%
      arrange(desc(key))

    # Classes with no segments anywhere are dropped from every panel rather
    # than drawn as an empty facet: with all values zero, ggplot pads the axis
    # to an arbitrary +/- range that reads like missing data instead of a
    # measured absence. The fact is stated in the subtitle instead.
    zero_cls <- comb %>% group_by(CLASS) %>%
      summarise(mx = max(FROH), .groups = "drop") %>% filter(mx == 0) %>% pull(CLASS)
    keep_cls <- setdiff(classes, zero_cls)
    keep_lab <- lab_cls[match(keep_cls, classes)]
    zero_note <- if (length(zero_cls))
      paste0(" | no segments detected in any individual for: ",
             paste(lab_cls[match(zero_cls, classes)], collapse = ", ")) else ""

    cb <- comb %>% filter(CLASS %in% keep_cls) %>%
      mutate(IID = factor(IID, levels = ord$IID),
             CLASS = factor(CLASS, levels = keep_cls, labels = keep_lab))

    p_res <- ggplot(cb, aes(IID, FROH, fill = DATASET)) +
      geom_col(position = position_dodge(preserve = "single"), width = 0.75) +
      facet_wrap(~ CLASS, ncol = 1, scales = "free_y") +
      scale_fill_manual(values = setNames(c("#f0a848", "#1a3a6b"), dcols_tag)) +
      labs(title = "F_ROH per individual, by dataset and ROH length class",
           subtitle = paste0("MAF >= 0.05 | ", dcols_tag[1], " = 9x target, ",
                             dcols_tag[2], " = 20x target | bars absent where the ",
                             "individual is not in that dataset",
                             sub("^ \\| ", "\n", zero_note)),
           x = "Individual (ordered by 20x short-class F_ROH)",
           y = expression(F[ROH]), fill = "Dataset") +
      theme_bw(base_size = 12) +
      theme(legend.position = "top", panel.grid.major.x = element_blank())
    save_plot(file.path(cdir, "froh_by_individual_both_datasets"), p_res, 11, 8)

    # ── PLOT 2: paired shift, shared individuals only ───────────────────────
    # A dumbbell makes the systematic direction visible in a way a scatter
    # does not: every segment for the short class points the same way.
    shared <- wide %>% filter(CLASS %in% keep_cls) %>%
      mutate(CLASS = factor(CLASS, levels = keep_cls, labels = keep_lab),
             IID = factor(IID, levels = rev(ord$IID[ord$IID %in% wide$IID])))
    p_dumb <- ggplot(shared, aes(y = IID)) +
      geom_segment(aes(x = .data[[dcols[1]]], xend = .data[[dcols[2]]],
                       yend = IID), colour = "grey60", linewidth = 0.8) +
      geom_point(aes(x = .data[[dcols[1]]], colour = dcols_tag[1]), size = 3) +
      geom_point(aes(x = .data[[dcols[2]]], colour = dcols_tag[2]), size = 3) +
      facet_wrap(~ CLASS, ncol = 3, scales = "free_x") +
      scale_colour_manual(values = setNames(c("#f0a848", "#1a3a6b"), dcols_tag)) +
      labs(title = "Change in F_ROH when the same individual is called at 20x instead of 9x",
           subtitle = paste0(nrow(shared) / length(keep_cls),
                             " individuals present in both datasets; every short-class segment points left",
                             sub("^ \\| ", "\n", zero_note)),
           x = expression(F[ROH]), y = NULL, colour = "Dataset") +
      theme_bw(base_size = 12) + theme(legend.position = "top")
    save_plot(file.path(cdir, "froh_dataset_shift"), p_dumb, 11, 5)

    # ── PLOT 3: concordance, with HONEST axes ───────────────────────────────
    # Equal x and y ranges per facet plus a square aspect ratio, so distance
    # from the dashed x=y line is read directly off the plot. With free
    # scales the axes stretch independently and points can look close to the
    # line when they are not.
    # Factored so facets read short -> medium -> long, not alphabetically.
    nz <- wide %>% filter(CLASS %in% keep_cls) %>%
      mutate(CLASS = factor(CLASS, levels = keep_cls, labels = keep_lab))
    if (length(zero_cls)) message("Plots: omitting all-zero class(es): ",
                                  paste(zero_cls, collapse = ", "))
    if (nrow(nz) > 0) {
      lims <- nz %>% group_by(CLASS) %>%
        summarise(m = max(c(.data[[dcols[1]]], .data[[dcols[2]]])) * 1.08, .groups = "drop")
      blank <- bind_rows(
        lims %>% mutate(x = 0, y = 0),
        lims %>% mutate(x = m, y = m))
      names(blank)[names(blank) == "x"] <- dcols[1]
      names(blank)[names(blank) == "y"] <- dcols[2]

      pc <- ggplot(nz, aes(.data[[dcols[1]]], .data[[dcols[2]]], label = IID)) +
        # inherit.aes = FALSE: the base mapping carries label = IID, which does
        # not exist in the limits frame and would error at build time.
        geom_blank(data = blank, inherit.aes = FALSE,
                   aes(x = .data[[dcols[1]]], y = .data[[dcols[2]]])) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
        geom_point(colour = "#1a3a6b", size = 3) +
        facet_wrap(~ CLASS, scales = "free") +
        labs(title = "F_ROH concordance between depth-equalised datasets",
             subtitle = paste0("Same individual called twice. Axes share a range within each panel, ",
                               "so distance from the\ndashed x=y line is read directly. ",
                               "Points below the line: the 9x call overstates F_ROH.",
                               sub("^ \\| ", "\n", zero_note)),
             x = paste0("F_ROH (", dcols_tag[1], ", 9x target)"),
             y = paste0("F_ROH (", dcols_tag[2], ", 20x target)")) +
        theme_bw(base_size = 12) +
        theme(legend.position = "none", aspect.ratio = 1)
      pc <- if (has_repel) pc + ggrepel::geom_text_repel(size = 3, max.overlaps = 20)
            else           pc + geom_text(size = 2.5, vjust = -0.5)
      save_plot(file.path(cdir, "froh_dataset_concordance"), pc, 10, 5.5)
    }
    message("Wrote: ", cdir)
  }
}
RSCRIPT
}

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────
declare -a TARGETS=()
if (( $# > 0 )); then
    for want in "$@"; do
        found=0
        for spec in "${DATASETS[@]}"; do
            [[ "${spec%%:*}" == "${want}" ]] && { TARGETS+=("${spec}"); found=1; }
        done
        (( found )) || die "unknown dataset '${want}' (known: ${DATASETS[*]%%:*})"
    done
else
    TARGETS=("${DATASETS[@]}")
fi

log "============================================================"
log " ROH analysis — ${#TARGETS[@]} dataset(s): ${TARGETS[*]%%:*}"
log " Genome size: ${GENOME_SIZE} bp"
log "============================================================"

preflight

for spec in "${TARGETS[@]}"; do
    tag="${spec%%:*}"; target="${spec##*:}"
    log "── ${tag} (target ${target}x) ──"
    density_check "${tag}"
    detect_roh "${tag}"
done

log "F_ROH calculation and plots"
R_SCRIPT="$(mktemp "${TMPDIR:-/tmp}/froh.XXXXXX.R")"
trap 'rm -f "${R_SCRIPT}"' EXIT
write_r_script "${R_SCRIPT}"

tag_list=(); for spec in "${TARGETS[@]}"; do tag_list+=("${spec%%:*}"); done
"${RSCRIPT_BIN}" "${R_SCRIPT}" "${OUT_DIR}" "${GENOME_SIZE}" "${PRIMARY_MAF_TAG}" "${tag_list[@]}" \
    2>&1 | tee -a "${MAIN_LOG}"

log "============================================================"
log " Complete"
for spec in "${TARGETS[@]}"; do
    tag="${spec%%:*}"
    log "   ${OUT_DIR}/${tag}/roh_analysis/froh_results/"
done
(( ${#TARGETS[@]} >= 2 )) && log "   ${OUT_DIR}/roh_comparison/    (cross-dataset concordance)"
log "============================================================"
