#!/bin/bash
# ---------------------------------------------------------------
# STEP 6 (pipeline v4). Runs of homozygosity on the FINAL files (mask v2)
# Based on Tatiana's ROH_tatiana.sh (21 Oct 2026); her logic and parameters
# are unchanged. Changes (marked "# J:"):
#   - inputs point at the local masked v2 all-sites (renamed on delivery)
#   - md5 and callable bp are hard checks (stop if they differ from hers)
#   - tool PATH set here (plink 1.9 before plink2)
#   - trees with no ROH in a class get an explicit 0 row
# Run from snp_calling/:  bash roh_v4.sh
# ---------------------------------------------------------------
set -euo pipefail

# J: tools. plink 1.9 + vcftools (popgen_tools) first, then bcftools/plink2
#    (SNP_call), bedtools (paralog_masks).
E=/home/genomics/miniconda3/envs
export PATH=$E/popgen_tools/bin:$E/SNP_call/bin:$E/paralog_masks/bin:$PATH
plink --version | grep -q 'v1.9' || { echo "plink is not 1.9" >&2; exit 1; }

# 0. Inputs: the masked v2 all-sites VCFs delivered on 1 Oct 2026
#    J: local names (delivered as callicarpa_d20_20x_9ind_allsites_masked_v2.vcf.gz
#    and callicarpa_d09_9x_13ind_allsites_masked_v2.vcf.gz)
A20=results/d20/allsites/cohort.allsites.masked_v2.vcf.gz       # md5 48435089038641b6a439c4ae00c7f238
A09=results/d09/allsites/cohort.allsites.masked_v2.vcf.gz       # md5 17bb7afe18f047c8d2bdb79b4a503da8
# J: hard check, both must match Tatiana's md5
# J: ONLY_FROH=1 reruns 6c-6d on existing outputs (skips md5, 6a, 6b)
ONLY_FROH=${ONLY_FROH:-0}
[ "$ONLY_FROH" = 1 ] || md5sum -c - <<EOF
48435089038641b6a439c4ae00c7f238  $A20
17bb7afe18f047c8d2bdb79b4a503da8  $A09
EOF
CHR=$(seq -f "Chr%02g" 1 17 | paste -sd,)                        # Chr01,...,Chr17 (scaffolds out)
declare -A EXPECTED=([d20]=159667585 [d09]=137293324)           # J: her callable bp

# ---------------------------------------------------------------
# 6a. ROH input: biallelic SNPs polymorphic among the trees, Chr01-Chr17,
#     90% call rate. No MAF and no HWE filter.
# ---------------------------------------------------------------
[ "$ONLY_FROH" = 1 ] || for D in d20 d09; do
  A=$([ $D = d20 ] && echo $A20 || echo $A09)
  mkdir -p v4/roh/$D
  # Polymorphic biallelic SNPs only (0 < AC < AN): the all-sites file also holds
  # ~3.5 M SNP records where every tree is homozygous ALT, which would lengthen ROH
  bcftools view -r $CHR -m2 -M2 -v snps -i 'INFO/AC>0 && INFO/AC<INFO/AN' \
      $A -Oz -o v4/roh/$D/snps.vcf.gz
  # 90% call rate, as in the old run
  vcftools --gzvcf v4/roh/$D/snps.vcf.gz --max-missing 0.9 \
      --recode --recode-INFO-all --stdout | bgzip > v4/roh/$D/roh_input.vcf.gz
  # J: SNP counts per step, for the Methods
  echo "$D SNPs: polymorphic biallelic $(bcftools view -H v4/roh/$D/snps.vcf.gz | wc -l)," \
       "after 90% call rate $(bcftools view -H v4/roh/$D/roh_input.vcf.gz | wc -l)" \
       | tee v4/roh/$D/snp_counts.txt
  # PLINK binary files (plink2 only converts; --homozyg needs PLINK 1.9)
  plink2 --vcf v4/roh/$D/roh_input.vcf.gz --allow-extra-chr \
      --set-missing-var-ids @:# --make-bed --out v4/roh/$D/roh
  # Callable sites as BED (every site of the all-sites VCF on Chr01-Chr17)
  bcftools query -r $CHR -f '%CHROM\t%POS0\t%POS\n' $A \
    | bedtools merge > v4/roh/$D/callable.bed
  L=$(awk '{s += $3 - $2} END {print s}' v4/roh/$D/callable.bed)
  echo "$D callable bp: $L (expected ${EXPECTED[$D]})"
  # J: hard check
  [ "$L" = "${EXPECTED[$D]}" ] || { echo "callable bp differs from Tatiana's" >&2; exit 1; }
done

# ---------------------------------------------------------------
# 6b. PLINK 1.9 --homozyg, three length classes (same as roh_claude.sh)
# ---------------------------------------------------------------
[ "$ONLY_FROH" = 1 ] || for D in d20 d09; do
  B=v4/roh/$D/roh
  #        class  min kb  min SNPs  max gap kb  window SNPs
  for p in "short 100 25 1000 25" "medium 500 50 1000 50" "long 1500 100 500 100"; do
    set -- $p
    plink --bfile $B --allow-extra-chr --homozyg \
      --homozyg-kb $2 --homozyg-snp $3 --homozyg-gap $4 --homozyg-density 50 \
      --homozyg-window-snp $5 --homozyg-window-het 1 --homozyg-window-missing 5 \
      --out v4/roh/$D/$1
  done
done

# ---------------------------------------------------------------
# 6c. Callable bp inside ROH, per tree and class -> F_ROH = callable bp in ROH / callable bp
# ---------------------------------------------------------------
for D in d20 d09; do
  L=$(awk '{s += $3 - $2} END {print s}' v4/roh/$D/callable.bed)
  echo "tree,class,n_seg,roh_callable_bp,froh" > v4/roh/$D/froh_per_tree.csv
  for k in short medium long; do
    H=v4/roh/$D/$k.hom
    # segments per tree (.hom: FID IID PHE CHR SNP1 SNP2 POS1 POS2 KB ...)
    awk 'NR > 1 {n[$2]++} END {for (t in n) print t"\t"n[t]}' $H > v4/roh/$D/$k.nseg
    # segments as BED, cut to callable sites, summed per tree
    # J: every tree of the .fam gets a row (0 when it has no ROH in this class)
    # J: plink2 recodes Chr01..Chr17 as 1..17; map back to the VCF names
    awk 'NR > 1 {c = ($4 ~ /^[0-9]+$/) ? sprintf("Chr%02d", $4) : $4; print c"\t"$7-1"\t"$8"\t"$2}' $H \
      | bedtools intersect -a - -b v4/roh/$D/callable.bed \
      | awk -v k=$k -v L=$L 'FILENAME == ARGV[1] {t[++nt] = $2; next}
             FILENAME == ARGV[2] {n[$1] = $2; next}
             {bp[$4] += $3 - $2}
             END {for (i = 1; i <= nt; i++) printf "%s,%s,%d,%d,%.6f\n", t[i], k, n[t[i]], bp[t[i]], bp[t[i]]/L}' \
            v4/roh/$D/roh.fam v4/roh/$D/$k.nseg - >> v4/roh/$D/froh_per_tree.csv
  done
done

# ---------------------------------------------------------------
# 6d. Versions for the Methods
# ---------------------------------------------------------------
{ bcftools --version | head -1; vcftools --version 2>&1 | head -1
  plink2 --version; plink --version; bedtools --version; } > v4/roh/versions.txt
# Send back: v4/roh/d20/ and v4/roh/d09/ (froh_per_tree.csv, short/medium/long .hom,
# .hom.indiv, .log, callable.bed) and v4/roh/versions.txt. snps.vcf.gz and roh_input.vcf.gz not needed.
