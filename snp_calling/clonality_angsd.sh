#!/bin/bash
# ---------------------------------------------------------------
# clonality_angsd.sh  (Part A)
# Genotype likelihoods (ANGSD) + pairwise IBD and relatedness (NgsRelate)
# on the native-depth BAMs of all 18 libraries. HEAVY: reads every BAM.
#
# Based on Tatiana Arias's script (Callicarpa_clonality_instructions_Jorge.txt,
# 5 Oct 2026). Her sites, read filters and ANGSD/NgsRelate options are kept.
# Changes, for robustness and speed:
#   * sites carry the called alleles (4 columns) and ANGSD uses them
#     (-doMajorMinor 3); the 18 positions where d09 and d20 disagree on the
#     alleles are dropped
#   * ANGSD runs in parallel by chromosome (17 chromosomes + 1 job for all
#     scaffolds); glf and mafs are concatenated in reference order
#   * the IBS matrix is rebuilt from the per-site .ibs.gz of all chunks, with
#     the number of shared sites per pair (per-chunk matrices can't be averaged)
#   * extra: NgsRelate rerun with allele frequencies estimated WITHOUT the two
#     repeat libraries 563 and 564 (her caveat about the frequencies)
#
# Usage (from /dados04/jorge/CALLICARPA_2/snp_calling, env "clonality"):
#   bash clonality_angsd.sh
#   TEST_REGION=Chr17:1-5000000 bash clonality_angsd.sh   # -> results/clonality_test
# ---------------------------------------------------------------
set -euo pipefail

# 0. Paths and resources
WORK_DIR=/dados04/jorge/CALLICARPA_2/snp_calling
REF=/dados04/jorge/CALLICARPA_2/CallicarpaGenome/car_asm.fa
TEST_REGION="${TEST_REGION:-}"
OUT=${WORK_DIR}/results/clonality${TEST_REGION:+_test}
JOBS=${JOBS:-18}                 # ANGSD jobs in parallel (one per chunk)
THREADS_PER_JOB=${THREADS_PER_JOB:-4}
THREADS=${THREADS:-16}           # NgsRelate
RSCRIPT=${RSCRIPT:-/home/genomics/miniconda3/envs/R_popstat_jorge/bin/Rscript}
D09=${WORK_DIR}/results/d09/masked/cohort.snps.masked_v2_gtfilt.vcf.gz
D20=${WORK_DIR}/results/d20/masked/cohort.snps.masked_v2_gtfilt.vcf.gz
cd "${WORK_DIR}"
for t in angsd ngsRelate bcftools parallel; do command -v $t >/dev/null || { echo "missing $t (conda activate clonality)" >&2; exit 1; }; done
mkdir -p $OUT/chunks
log() { echo "[$(date +'%F %T')] $*"; }

# 1. Sites = union of the masked SNP positions in d09 and d20, with alleles,
#    in reference (.fai) order. The paralog mask and genotype filters apply.
log "sites"
( bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $D09
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $D20 ) \
  | awk -F'\t' -v OFS='\t' -v conf=$OUT/sites_allele_conflicts.txt '
      NR==FNR { o[$1]=NR; next }
      { k=$1":"$2; a=$3"\t"$4
        if (!(k in al)) { al[k]=a; ord[k]=o[$1]"\t"$1"\t"$2 }
        else if (al[k]!=a) bad[k]=1 }
      END { for (k in al) if (!(k in bad)) print ord[k], al[k]
            for (k in bad) print k > conf }' $REF.fai - \
  | sort -k1,1n -k3,3n | cut -f2- > $OUT/sites.all.txt
touch $OUT/sites_allele_conflicts.txt
if [[ -n "$TEST_REGION" ]]; then
  c=${TEST_REGION%%:*}; r=${TEST_REGION#*:}; s=${r%-*}; e=${r#*-}
  awk -v c=$c -v s=$s -v e=$e '$1==c && $2>=s && $2<=e' $OUT/sites.all.txt > $OUT/sites.txt
else
  mv $OUT/sites.all.txt $OUT/sites.txt
fi
rm -f $OUT/sites.all.txt
# ANGSD refuses an index that is not strictly newer than sites.txt: rebuild it
rm -f $OUT/sites.txt.bin $OUT/sites.txt.idx; sleep 2
angsd sites index $OUT/sites.txt 2> $OUT/sites_index.log
log "sites: $(wc -l < $OUT/sites.txt) (dropped for d09/d20 allele conflict: $(wc -l < $OUT/sites_allele_conflicts.txt))"

# 2. BAM list in a fixed order (547 ... 564) and matching ids for NgsRelate
ls results/native/bam/*.markdup.bam | sort > $OUT/bams.txt
sed 's#.*/##; s#\.markdup\.bam##' $OUT/bams.txt > $OUT/ids.txt
[ $(wc -l < $OUT/bams.txt) -eq 18 ] || { echo "expected 18 BAMs" >&2; exit 1; }
N=18

# 3. Genotype likelihoods + per-site IBS (ANGSD), one job per chunk.
#    Read filters (Tatiana): unique hits, proper pairs, MAPQ >= 30, base quality
#      >= 20, BAQ, -remove_bads drops duplicate-flagged reads
#    -GL 2 GATK model; -doGlf 3 binary likelihoods for NgsRelate
#    -doMajorMinor 3: alleles from the sites file; -doMaf 1: knownEM frequencies
#    -SNP_pval 1e-6 -minMaf 0.05: variable in these 18 libraries
#    -doIBS 1 -output01 1: one random read per sample per site, written per site
#    (Frequency-free genotype calls are made later in R from the same glf,
#     not with -doPost/-doGeno: those options change which sites ANGSD keeps.)
#    -minInd 9: data in at least 9 of the 18 libraries
if [[ -n "$TEST_REGION" ]]; then
  printf 'test\t-r\t%s\n' "$TEST_REGION" > $OUT/chunks.tsv
else
  cut -f1 $OUT/sites.txt | uniq > $OUT/regions.txt          # contigs that carry sites
  grep '^Chr' $OUT/regions.txt | awk -v OFS='\t' '{print $1, "-r", $1":"}' > $OUT/chunks.tsv
  grep -v '^Chr' $OUT/regions.txt > $OUT/chunks/scaffolds.rf || true
  [ -s $OUT/chunks/scaffolds.rf ] && printf 'scaffolds\t-rf\t%s\n' "$OUT/chunks/scaffolds.rf" >> $OUT/chunks.tsv
fi
run_chunk() {
  set -euo pipefail
  local name=$1 flag=$2 reg=$3 seed=$4
  if [[ -f $OUT/chunks/$name.done ]]; then echo "chunk $name: done, skipping"; return 0; fi
  angsd -b $OUT/bams.txt -ref $REF -out $OUT/chunks/$name \
    -sites $OUT/sites.txt $flag $reg \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
    -minMapQ 30 -minQ 20 -C 50 -baq 1 \
    -GL 2 -doGlf 3 -doMajorMinor 3 -doMaf 1 -SNP_pval 1e-6 -minMaf 0.05 \
    -doCounts 1 -doIBS 1 -output01 1 \
    -minInd 9 -seed $seed -P $THREADS_PER_JOB > $OUT/chunks/$name.angsd.log 2>&1
  touch $OUT/chunks/$name.done
  echo "[$(date +'%F %T')] chunk $name: $(zcat $OUT/chunks/$name.mafs.gz | sed 1d | wc -l) sites"
}
export -f run_chunk; export OUT REF THREADS_PER_JOB
log "ANGSD: $(wc -l < $OUT/chunks.tsv) chunks, $JOBS in parallel x $THREADS_PER_JOB threads"
awk -v OFS='\t' '{print $0, 42+NR}' $OUT/chunks.tsv \
  | PARALLEL_SHELL=/bin/bash parallel --colsep '\t' --halt now,fail=1 -j $JOBS run_chunk {1} {2} {3} {4}

# Merge in reference order: glf (gzip streams concatenate), mafs (one header), ibs
names=$(cut -f1 $OUT/chunks.tsv)
cat $(for n in $names; do echo $OUT/chunks/$n.glf.gz; done) > $OUT/all18.glf.gz
first=$(awk 'NR==1{print $1}' $OUT/chunks.tsv)
{ zcat $OUT/chunks/$first.mafs.gz | awk 'NR==1'
  for n in $names; do zcat $OUT/chunks/$n.mafs.gz | sed 1d; done; } | gzip > $OUT/all18.mafs.gz
{ zcat $OUT/chunks/$first.ibs.gz | awk 'NR==1'
  for n in $names; do zcat $OUT/chunks/$n.ibs.gz | sed 1d; done; } \
  | sed 's/\t$//' | gzip > $OUT/all18.ibs.gz           # ANGSD ends data rows with a tab
nsites=$(zcat $OUT/all18.mafs.gz | sed 1d | wc -l)
glfbytes=$(zcat $OUT/all18.glf.gz | wc -c)
[ "$glfbytes" -eq $((nsites * N * 3 * 8)) ] || { echo "glf size $glfbytes != $nsites sites x $N x 3 doubles" >&2; exit 1; }
{ cat $OUT/chunks/$first.arg
  echo; echo "# Parallel run by chunk (clonality_angsd.sh). Sites used per chunk:"
  for n in $names; do echo "$n	$(zcat $OUT/chunks/$n.mafs.gz | sed 1d | wc -l)"; done
  echo "TOTAL	$nsites"; } > $OUT/all18.arg
log "merged: $nsites sites in all18.glf.gz / all18.mafs.gz"

# 4. NgsRelate: J9, J8, J7 and rab for every pair.
#    Frequencies = knownEM column of the .mafs file (found by its header name)
zcat $OUT/all18.mafs.gz | awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="knownEM") c=i; if (!c) exit 1; next} {print $c}' > $OUT/freq.txt
ngsRelate -g $OUT/all18.glf.gz -n $N -f $OUT/freq.txt -z $OUT/ids.txt -O $OUT/ngsrelate.res -p $THREADS > $OUT/ngsrelate.log 2>&1
log "NgsRelate: $(($(wc -l < $OUT/ngsrelate.res) - 1)) pairs"

# 5. R: IBS matrix + shared sites, UPGMA png, frequencies without 563/564
#    (glf is BGZF: read it through zcat; R's gzcon stops after the first block)
"$RSCRIPT" - "$OUT" $N <<'REOF'
a <- commandArgs(TRUE); out <- a[1]; n <- as.integer(a[2])
ids <- readLines(file.path(out, "ids.txt"))
zc <- function(f) pipe(paste("zcat", shQuote(file.path(out, f))))
# IBS from per-site random reads (-output01 1: 0/1 allele, -1 missing)
ibs <- read.delim(zc("all18.ibs.gz"))
g <- as.matrix(ibs[, grep("^ind", names(ibs))]); g[g < 0] <- NA
D <- matrix(0, n, n, dimnames = list(ids, ids)); S <- D
for (i in 1:(n - 1)) for (j in (i + 1):n) {
  ok <- !is.na(g[, i]) & !is.na(g[, j])
  S[i, j] <- S[j, i] <- sum(ok)
  D[i, j] <- D[j, i] <- mean(g[ok, i] != g[ok, j])
}
write.table(D, file.path(out, "all18.ibsMat"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
png(file.path(out, "upgma_ibs.png"), width = 1000, height = 650, res = 120)
plot(hclust(as.dist(D), method = "average"), main = "UPGMA of IBS distance (18 libraries)",
     xlab = "", sub = "", ylab = "IBS distance"); invisible(dev.off())
# Allele frequencies by EM from the glf, without the repeat libraries 563/564
ns <- length(readLines(file.path(out, "freq.txt")))
con <- zc("all18.glf.gz"); open(con, "rb"); gl <- readBin(con, "double", n = ns * n * 3 + 1); close(con)
stopifnot(length(gl) == ns * n * 3)
L <- exp(array(gl, dim = c(3, n, ns)))            # log-likelihoods -> likelihoods
# Frequency-free check: genotypes called with a FLAT prior (= ANGSD -doPost 2),
# kept if the posterior is >= 0.95; libraries without reads stay missing
Sum <- L[1, , ] + L[2, , ] + L[3, , ]
P0 <- L[1, , ] / Sum; P1 <- L[2, , ] / Sum; P2 <- L[3, , ] / Sum
best <- ifelse(P2 >= P1 & P2 >= P0, 2, ifelse(P1 >= P0, 1, 0))
best[pmax(P0, P1, P2) < 0.95] <- NA
G <- t(best)                                       # sites x libraries
em <- function(keep) {
  k <- length(keep); Lk <- L[, keep, , drop = FALSE]; p <- rep(0.25, ns)
  for (it in 1:500) {
    w0 <- (1 - p)^2; w1 <- 2 * p * (1 - p); w2 <- p^2
    num <- Lk[2, , ] * rep(w1, each = k) + 2 * Lk[3, , ] * rep(w2, each = k)
    den <- Lk[1, , ] * rep(w0, each = k) + Lk[2, , ] * rep(w1, each = k) + Lk[3, , ] * rep(w2, each = k)
    pn <- colSums(num / den) / (2 * k)
    if (max(abs(pn - p)) < 1e-7) { p <- pn; break }
    p <- pn
  }
  p
}
f18 <- scan(file.path(out, "freq.txt"), quiet = TRUE); p18 <- em(seq_len(n))
cat(sprintf("EM check, 18 libraries vs ANGSD knownEM: max |diff| = %.2g, cor = %.6f\n", max(abs(p18 - f18)), cor(p18, f18)))
p16 <- em(which(!ids %in% c("563", "564")))
writeLines(format(p16, digits = 6, scientific = FALSE), file.path(out, "freq16.txt"))
# Per-pair table: IBS + called-genotype discordance + opposite homozygotes
pr <- t(combn(n, 2))
tab <- data.frame(ida = ids[pr[, 1]], idb = ids[pr[, 2]], ibs_sites = S[pr], ibs_dist = round(D[pr], 5))
tab$geno_sites <- tab$geno_discord <- tab$opp_hom <- NA
for (k in seq_len(nrow(pr))) {
  x <- G[, pr[k, 1]]; y <- G[, pr[k, 2]]; ok <- !is.na(x) & !is.na(y)
  tab$geno_sites[k] <- sum(ok)
  tab$geno_discord[k] <- round(mean(x[ok] != y[ok]), 5)
  tab$opp_hom[k] <- round(mean(abs(x[ok] - y[ok]) == 2), 5)
}
write.table(tab, file.path(out, "pairs_ibs_geno.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
REOF
ngsRelate -g $OUT/all18.glf.gz -n $N -f $OUT/freq16.txt -z $OUT/ids.txt -O $OUT/ngsrelate_freq16.res -p $THREADS > $OUT/ngsrelate_freq16.log 2>&1

# 6. One table per pair: NgsRelate (both frequency sets) + IBS + called genotypes.
#    Classes follow Tatiana's table; the repeat-library pairs are listed both
#    ways (562-563/561-564 as in her table, 561-563/562-564 as the data suggest).
"$RSCRIPT" - "$OUT" <<'REOF'
out <- commandArgs(TRUE)[1]
r18 <- read.delim(file.path(out, "ngsrelate.res")); r16 <- read.delim(file.path(out, "ngsrelate_freq16.res"))
x <- read.delim(file.path(out, "pairs_ibs_geno.tsv"))
key <- function(d) paste(d$ida, d$idb)
k18 <- r18[, c("ida", "idb", "nSites", "J9", "J8", "J7", "rab")]
k16 <- r16[, c("ida", "idb", "J7", "rab")]; names(k16)[3:4] <- c("J7_no563_564", "rab_no563_564")
m <- merge(merge(k18, k16, by = c("ida", "idb")), x, by = c("ida", "idb"))
a <- m$ida; b <- m$idb
m$class <- ifelse(paste(a, b) %in% c("562 563", "561 564"), "control_as_listed",
           ifelse(paste(a, b) %in% c("561 563", "562 564"), "control_swapped",
           ifelse(a == 549 & b %in% 550:552, "549_vs_PV",
           ifelse(a == 549 & b %in% 553:560, "549_vs_Sonadora", "other"))))
m <- m[order(m$class != "control_as_listed", m$class != "control_swapped", m$class != "549_vs_PV",
             m$class != "549_vs_Sonadora", m$geno_discord), ]
write.table(m, file.path(out, "pairs_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
poi <- rbind(m[m$class != "other", ], head(m[m$class == "other", ], 5))
write.table(poi, file.path(out, "pairs_of_interest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("== controls (same tree expected: J7, rab ~ 1; genotype discordance near the error rate)\n")
print(m[grepl("^control", m$class), c("ida", "idb", "nSites", "J7", "rab", "rab_no563_564", "ibs_dist", "geno_sites", "geno_discord", "opp_hom")], row.names = FALSE)
REOF

# 7. Versions for the Methods
{ angsd 2>&1 | awk '/angsd version/ && !f {print; f=1}'
  ngsRelate 2>&1 | awk '/NgsRelate version/ && !f {print; f=1}'
  bcftools --version | awk 'NR==1'; } > $OUT/versions.txt
log "done -> $OUT"
