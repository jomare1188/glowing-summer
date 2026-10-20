# Relative-depth histograms for the ExcessHet layer (mask_v3_genets.sh reldepth).
# Input dir holds reldepth_bins.tsv (0.01-wide bins, 6 = ">= 6x") and
# reldepth_groups.tsv. Writes reldepth_summary.tsv and reldepth_hist.png/.pdf.
# A collapsed 2-copy paralog is expected near 2x, a real het near 1x.

args <- commandArgs(trailingOnly = TRUE)
d <- args[1]
b <- read.delim(file.path(d, "reldepth_bins.tsv"))
g <- read.delim(file.path(d, "reldepth_groups.tsv"))

groups <- c("not_flagged", "flagged", "flagged_bal_allhet", "flagged_rest")
groups <- groups[groups %in% b$group]

# ---- summary from the 0.01 bins (median / mode at 0.01x resolution) ----
summ <- function(x) {
  x <- x[order(x$bin_lo), ]
  n <- sum(x$n); cs <- cumsum(x$n)
  med <- x$bin_lo[which(cs >= n / 2)[1]] + 0.005
  # mode on 0.05 bins inside 0.1-4x (skips the zero-depth spike)
  b5 <- tapply(x$n, floor(x$bin_lo / 0.05 + 1e-9) * 0.05, sum)
  b5 <- b5[as.numeric(names(b5)) >= 0.1 & as.numeric(names(b5)) < 4]
  mode <- as.numeric(names(b5)[which.max(b5)]) + 0.025
  p <- function(lo, hi) round(100 * sum(x$n[x$bin_lo >= lo - 1e-9 & x$bin_lo < hi - 1e-9]) / n, 1)
  c(n_sites = n, median_x = round(med, 3), mode_x = round(mode, 3),
    pct_lt0.5 = p(0, 0.5), pct_0.5_1.5 = p(0.5, 1.5), pct_1.5_2.5 = p(1.5, 2.5), pct_ge2.5 = p(2.5, Inf))
}
rows <- list()
for (nm in c("cutoffs", "vcf")) for (gr in groups) {
  s <- summ(b[b$group == gr & b$norm == nm, ])
  hd <- g$n_in_v2_highdepth[g$group == gr]
  rows[[length(rows) + 1]] <- data.frame(normalisation = nm, group = gr, t(s),
                                         pct_in_v2_highdepth = round(100 * hd / s[["n_sites"]], 1),
                                         check.names = FALSE)
}
out <- do.call(rbind, rows)
write.table(out, file.path(d, "reldepth_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- histograms: proportion of sites per 0.05x bin, 0-4x ----
dens <- function(gr, nm) {
  x <- b[b$group == gr & b$norm == nm, ]
  tot <- sum(x$n)
  x <- x[x$bin_lo < 4, ]
  h <- tapply(x$n, floor(x$bin_lo / 0.05 + 1e-9) * 0.05, sum)
  br <- seq(0, 3.95, by = 0.05)
  y <- setNames(rep(0, length(br)), sprintf("%.2f", br))
  y[sprintf("%.2f", as.numeric(names(h)))] <- h
  list(x = br, y = 100 * y / tot)
}
col <- c(not_flagged = "grey45", flagged = "firebrick", flagged_bal_allhet = "darkorange3")
lab <- function(gr, nm) {
  r <- out[out$group == gr & out$normalisation == nm, ]
  sprintf("%s (n = %s; median %.2fx)", gr, format(r$n_sites, big.mark = ","), r$median_x)
}
panel <- function(gs, nm, main) {
  ds <- lapply(gs, dens, nm = nm)
  ymax <- max(sapply(ds, function(z) max(z$y)))
  plot(NA, xlim = c(0, 4), ylim = c(0, ymax * 1.08), xaxs = "i", yaxs = "i",
       xlab = "Mean depth relative to each genet's median (7 genets)",
       ylab = "% of sites per 0.05x bin", main = main)
  abline(v = c(1, 2), lty = 2, col = "grey30")
  for (i in seq_along(gs)) lines(ds[[i]]$x + 0.025, ds[[i]]$y, type = "s", lwd = 2, col = col[gs[i]])
  legend("topright", legend = sapply(gs, lab, nm = nm), col = col[gs], lwd = 2, bty = "n", cex = 0.8)
}
draw <- function() {
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
  panel(c("not_flagged", "flagged"), "cutoffs", "A. Flagged (ExcessHet >= 20) vs not flagged\nmedians: cutoffs.tsv (mosdepth)")
  panel("flagged_bal_allhet", "cutoffs", "B. Flagged, balanced AB and all genets het\nmedians: cutoffs.tsv (mosdepth)")
  panel(c("not_flagged", "flagged"), "vcf", "A'. As A, medians: VCF FORMAT/DP\n(not-flagged SNPs; scale check)")
  panel("flagged_bal_allhet", "vcf", "B'. As B, medians: VCF FORMAT/DP\n(not-flagged SNPs; scale check)")
}
png(file.path(d, "reldepth_hist.png"), width = 2400, height = 1800, res = 200); draw(); invisible(dev.off())
pdf(file.path(d, "reldepth_hist.pdf"), width = 12, height = 9); draw(); invisible(dev.off())
