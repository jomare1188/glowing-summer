#!/bin/bash
################################################################################
# HETEROZYGOSITY CHECKS  (feedback.txt, Tatiana 28 Sep 2026)
#
# After the paralog masks, d20 still shows Fis -0.16 (about -0.06 is expected
# from n = 9 alone). Tatiana's checks, run identically for any VCF x mask:
#   1. het-count spectrum: SNPs with k = 0..n heterozygous trees (all trees
#      called), observed vs the exact HWE expectation given each site's allele count.
#      An excess at k = n-1, n means paralogs remain.
#   2. per tree: Ho, hom-alt rate, F (vcftools --het), share of skewed het
#      calls (alt fraction < 0.25 or > 0.75, depth >= 8) and the allele-balance
#      histogram. One or two outlier trees would mean a sample problem.
#   3. Fis: mean per-site Fis (as fis_before_after.tsv) and ratio-of-sums Fis
#      with He corrected for sample size (2n/(2n-1)), which is ~0 under HWE.
#
# Usage:  VCF=<vcf.gz> MASK=<exclusion bed|none> TAG=<name> \
#         [SAMPLES=<file>] [PREFILTER=hard] bash het_checks.sh
#   PREFILTER=hard applies the GATK hard filters of snp_datasets.sh with
#   bcftools, for the raw native-depth calls. Outputs: results/het_checks/.
#   PLOTS_ONLY=1 TAG=<name> bash het_checks.sh  redraws the plots of an
#   existing TAG from its tables, without recomputing.
################################################################################

set -euo pipefail

PROJECT_DIR="/dados04/jorge/CALLICARPA_2"
WORK_DIR="${PROJECT_DIR}/snp_calling"
OUT="${WORK_DIR}/results/het_checks"
TMP_DIR="${WORK_DIR}/tmp"
RSCRIPT_BIN="${RSCRIPT_BIN:-/home/genomics/miniconda3/envs/R_popstat_jorge/bin/Rscript}"
VCFTOOLS_BIN="${VCFTOOLS_BIN:-/home/genomics/miniconda3/envs/popgen_tools/bin/vcftools}"

PLOTS_ONLY="${PLOTS_ONLY:-0}"
VCF="${VCF:-}"
MASK="${MASK:-none}"
TAG="${TAG:?set TAG}"
SAMPLES="${SAMPLES:-}"
PREFILTER="${PREFILTER:-}"

# Same hard filters as HARD_FILTERS in snp_datasets.sh
HARD_EXPR='QD<2.0 || FS>60.0 || MQ<40.0 || SOR>3.0 || MQRankSum<-12.5 || ReadPosRankSum<-8.0'

work="${TMP_DIR}/het_checks_${TAG}"
if [[ "${PLOTS_ONLY}" != 1 ]]; then
[[ -n "${VCF}" ]] || { echo "set VCF" >&2; exit 1; }
command -v bcftools >/dev/null || { echo "conda activate SNP_call" >&2; exit 1; }
mkdir -p "${OUT}" "${TMP_DIR}" "${work}"
echo "[$(date +'%F %T')] ${TAG}: $(basename "${VCF}"), mask $(basename "${MASK}")${SAMPLES:+, samples $(basename "${SAMPLES}")}${PREFILTER:+, prefilter ${PREFILTER}}"

# -- 1. the SNP set: biallelic SNPs, masked, optional sample subset ------------
sel="${work}/snps.bcf"
{
    args=(-m2 -M2 -v snps)
    [[ "${MASK}" != none ]] && args+=(-T "^${MASK}")
    [[ -n "${SAMPLES}" ]] && args+=(-S "${SAMPLES}")
    bcftools view "${args[@]}" -Ou "${VCF}" \
      | { if [[ "${PREFILTER}" == hard ]]; then bcftools filter -e "${HARD_EXPR}" -Ou; else cat; fi; } \
      | { if [[ "${PREFILTER}" == hard ]]; then bcftools view -i 'F_MISSING<=0.2' -Ou; else cat; fi; } \
      | bcftools view -e 'INFO/AC==0 || INFO/AC==INFO/AN' -Ob -o "${sel}"
}
bcftools index -f "${sel}"
bcftools query -l "${sel}" > "${work}/samples.txt"
n=$(wc -l < "${work}/samples.txt")
echo "  $(bcftools index -n "${sel}") SNPs, ${n} trees"

# -- 2. one pass over genotypes + AD -------------------------------------------
bcftools query -f '[%GT:%AD\t]\n' "${sel}" \
  | awk -F'\t' -v OFS='\t' -v n="${n}" -v names="$(paste -sd, "${work}/samples.txt")" \
        -v spec="${OUT}/het_spectrum_${TAG}.tsv" -v ps="${work}/per_sample_core.tsv" \
        -v ab="${OUT}/allele_balance_${TAG}.tsv" -v fis="${work}/fis.tsv" '
    BEGIN{ split(names, nm, ","); lg[0]=0; for(i=1;i<=n;i++) lg[i]=lg[i-1]+log(i) }
    { nc=0; nh=0; alt=0
      for(i=1;i<=n;i++){ split($i,x,":"); g=x[1]; if(g ~ /\./) continue
        nc++; split(g,a,/[\/|]/); alt+=a[1]+a[2]; c[i]++
        if(a[1]!=a[2]){ nh++; het[i]++
          split(x[2],d,","); dp=d[1]+d[2]
          if(dp>=8){ f=d[2]/dp; nab[i]++; if(f<0.25||f>0.75) skew[i]++; b=int(f*20); if(b>19)b=19; hist[i,b]++ } }
        else if(a[1]==1) homalt[i]++ }
      if(nc<2) next
      p=alt/(2*nc); he=2*p*(1-p); if(he==0) next
      ho=nh/nc; ns++; sfis+=1-ho/he; sho+=ho; she+=he; shec+=he*2*nc/(2*nc-1)
      if(nc==n){ obs[nh]++; full++; nm_alt[alt]++ } }
    END{
      # Expected spectrum: exact HWE distribution of the het count GIVEN each
      # site alt-allele count m (Levene 1949; the distribution behind the exact
      # HWE test). Conditioning on m removes the plug-in bias of Binom(n, 2pq).
      for(i=n+1;i<=2*n;i++) lg[i]=lg[i-1]+log(i)
      for(m in nm_alt){ m+=0
        for(h=m%2; h<=m && h<=2*n-m; h+=2){ aa=(m-h)/2; rr=n-h-aa; if(rr<0) continue
          lp=lg[n]-lg[rr]-lg[h]-lg[aa]+h*log(2)+lg[m]+lg[2*n-m]-lg[2*n]
          ex[h]+=nm_alt[m]*exp(lp) } }
      print "k_het","observed","expected_hwe","obs_over_exp" > spec
      for(k=0;k<=n;k++) printf "%d\t%d\t%.1f\t%.3f\n", k, obs[k]+0, ex[k], (ex[k]>0?(obs[k]+0)/ex[k]:0) > spec
      print "sample","called","Ho","homalt_rate","het_calls_ab","frac_skewed_ab","median_alt_frac" > ps
      for(i=1;i<=n;i++){ # median alt fraction, from the 20-bin histogram (bin centre)
        cnt=0; for(b=0;b<20;b++){ cnt+=hist[i,b]; if(cnt>=nab[i]/2){ med=(b+0.5)/20; break } }
        printf "%s\t%d\t%.4f\t%.4f\t%d\t%.4f\t%.3f\n", nm[i], c[i], het[i]/c[i], homalt[i]/c[i], nab[i], skew[i]/nab[i], med > ps }
      print "sample","bin_low","frac" > ab
      for(i=1;i<=n;i++) for(b=0;b<20;b++) printf "%s\t%.2f\t%.5f\n", nm[i], b/20, hist[i,b]/nab[i] > ab
      printf "%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n", ns, full, sfis/ns, sho/ns, 1-sho/she, 1-sho/shec > fis }'

# -- 3. F per tree (vcftools --het) and outlier flag (> 3 MAD from median Ho) ---
bcftools view -Oz -o "${work}/snps.vcf.gz" "${sel}"
"${VCFTOOLS_BIN}" --gzvcf "${work}/snps.vcf.gz" --het --out "${work}/vt" 2> "${work}/vcftools.log"
awk -F'\t' -v OFS='\t' 'NR==FNR{ if(FNR>1) F[$1]=$5; next }
    FNR==1{ print $0, "F_vcftools"; next } { print $0, F[$1] }' \
    "${work}/vt.het" "${work}/per_sample_core.tsv" > "${work}/ps_F.tsv"
"${RSCRIPT_BIN}" -e '
  x <- read.delim(commandArgs(TRUE)[1]); m <- median(x$Ho); d <- mad(x$Ho)
  x$outlier_Ho_3MAD <- ifelse(abs(x$Ho - m) > 3 * d, "YES", "no")
  write.table(x, commandArgs(TRUE)[2], sep = "\t", quote = FALSE, row.names = FALSE)' \
  "${work}/ps_F.tsv" "${OUT}/per_sample_${TAG}.tsv"

# -- 4. Fis summary row (one row per TAG, replaced on rerun) --------------------
fsum="${OUT}/fis_summary.tsv"
[[ -f "${fsum}" ]] || printf 'tag\tn_trees\tsnps\tsnps_all_called\tmean_site_Fis\tmean_Ho\tFis_ratio\tFis_ratio_corrected\texpected_site_Fis\n' > "${fsum}"
grep -v -P "^${TAG}\t" "${fsum}" > "${fsum}.tmp" || true
awk -F'\t' -v OFS='\t' -v t="${TAG}" -v n="${n}" '{ print t, n, $1, $2, $3, $4, $5, $6, sprintf("%.4f", -1/(2*n-1)) }' "${work}/fis.tsv" >> "${fsum}.tmp"
mv "${fsum}.tmp" "${fsum}"

fi   # PLOTS_ONLY

# -- 5. plots -------------------------------------------------------------------
"${RSCRIPT_BIN}" - "${OUT}" "${TAG}" <<'REOF'
a <- commandArgs(TRUE); out <- a[1]; tag <- a[2]
sp <- read.delim(file.path(out, paste0("het_spectrum_", tag, ".tsv")))
ps <- read.delim(file.path(out, paste0("per_sample_", tag, ".tsv")))
ab <- read.delim(file.path(out, paste0("allele_balance_", tag, ".tsv")))
n <- max(sp$k_het)
dev <- function(f, w, h, expr) {
  png(file.path(out, paste0(f, "_", tag, ".png")), width = w * 120, height = h * 120, res = 120); expr(); invisible(dev.off())
  pdf(file.path(out, paste0(f, "_", tag, ".pdf")), width = w, height = h); expr(); invisible(dev.off()) }
dev("het_spectrum", 11, 4.5, function() {
  par(mfrow = c(1, 2), mar = c(4.5, 6.5, 3, 1))
  m <- rbind(sp$observed, sp$expected_hwe); colnames(m) <- sp$k_het
  barplot(m, beside = TRUE, col = c("grey25", "grey75"), border = NA, las = 1,
          xlab = paste0("heterozygous trees per SNP (all ", n, " called)"),
          main = paste0(tag, ": observed vs HWE"), legend.text = c("observed", "expected (HWE)"),
          args.legend = list(x = "topright", bty = "n"))
  title(ylab = "SNPs", line = 5)   # clear of the horizontal tick labels
  plot(sp$k_het, sp$obs_over_exp, type = "b", pch = 19, log = "y", las = 1,
       xlab = "heterozygous trees per SNP", ylab = "observed / expected", main = "excess by class")
  abline(h = 1, lty = 2, col = "firebrick") })
dev("per_sample", 8, 4.5, function() {
  par(mar = c(4.5, 5, 3, 1)); o <- order(ps$Ho)
  cols <- ifelse(ps$outlier_Ho_3MAD[o] == "YES", "firebrick", "grey30")
  plot(seq_along(o), ps$Ho[o], pch = 19, col = cols, xaxt = "n", las = 1, ylim = range(0, ps$Ho),
       xlab = "", ylab = "Ho (het / called SNPs)", main = paste0(tag, ": Ho per tree (red = >3 MAD)"))
  axis(1, at = seq_along(o), labels = ps$sample[o], las = 2); abline(h = median(ps$Ho), lty = 2)
  points(seq_along(o), ps$frac_skewed_ab[o], pch = 4, col = "steelblue")
  legend("bottomright", c("Ho", "share of skewed het calls"), pch = c(19, 4), col = c("grey30", "steelblue"), bty = "n") })
s <- unique(ab$sample); nc <- ceiling(length(s) / 3)
dev("allele_balance", 3.2 * nc, 8, function() {
  par(mfrow = c(3, nc), mar = c(3.5, 3.5, 2.5, 0.5))
  for (i in s) { z <- ab[ab$sample == i, ]
    barplot(z$frac, space = 0, col = "grey55", border = NA, main = i, las = 1, ylim = c(0, max(ab$frac)))
    axis(1, at = c(0, 5, 10, 15, 20), labels = c(0, .25, .5, .75, 1)); abline(v = 10, lty = 2) } })
REOF

rm -rf "${work}"
echo "[$(date +'%F %T')] ${TAG}: done -> ${OUT}/*_${TAG}.{tsv,png,pdf}, fis_summary.tsv"
