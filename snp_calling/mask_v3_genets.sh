#!/bin/bash
################################################################################
# PARALOG MASK v3: ONE LIBRARY PER GENET
#
# Why this exists
#   Clonality Part A (5 Oct 2026) showed that the 13 libraries used to build
#   the v2 mask (paralog_masks.sh) are only 7 genets:
#       Naguabo genet   547, 548, 549
#       FS831 genet     554, 555, 556
#       FS836s6 genet   559, 560, library 562
#       550 (PV2FS), 551 (PV3FS), 552 (PV4FS), library 561
#   Two layers of the v2 mask count trees, so they counted these genets 2-3x:
#     - ExcessHet >= 20 (mask_paralog_windows.bed, 21.2 % of the genome).
#       ExcessHet's ceiling with 7 independent genets is below 19.8 (the n = 9
#       value), so with genets counted once no site could reach 20. Every
#       flagged site passed only because clonal ramets added extra het calls.
#       ExcessHet cannot be used at n = 7; this script replaces it with read
#       allele balance plus an "all genets heterozygous" rule.
#     - High depth ">= 7 of 13 trees" (mask_highdepth.bed). Rebuilt here as
#       ">= 4 of 7 genets" from the per-tree BEDs that already exist.
#   The 1-kb window was chosen to bring Fis to the HWE value for 9 independent
#   trees; the d20 trees are 7 genets, so the window is re-chosen here.
#   Mappability (mask_lowmap.bed) depends only on the reference: unchanged.
#   The genotype filters (DP < 3, het binomial balance) are not part of the
#   mask: unchanged.
#
# One library per genet: the deepest ramet of each genet (depth_summary.tsv):
#   547 (41.7x), 554 (50.2x), 559 (67.7x), 550, 551, 552, 561.
#   All seven are also in d20, so the window test runs on d20 with them.
#
# Phases (run in this order; stop after "check" and agree thresholds):
#   stats    per-site H and read allele balance in the 7 genets, from the
#            native-depth joint calls; keeps the 13-library ExcessHet
#   check    QUICK CHECK for Tatiana: of the sites flagged by ExcessHet >= 20,
#            how many show paralog evidence in the 7 genets, and how many look
#            like real heterozygosity (balanced reads, not all genets het)
#   reldepth depth relative to each genet's median, flagged vs not flagged
#            (Tatiana 19 Oct: near 2x -> mask stays; near 1x -> go on)
#   ---- review the check with Tatiana; if the mask stays, stop here ----
#   depth    high-depth layer counting genets (>= HIGHDP_MIN_GENETS of 7)
#   paralog  new paralog sites: skewed allele balance OR all genets het
#   window   compares windows 0 / 500 / 1000 bp on d20 (7 genets) with
#            het_checks.sh; main criterion = share of skewed het calls
#   union    v3 union (lowmap + depth by genet + paralog windows) and how much
#            differs from v2
#   remask   masked all-sites and SNP VCFs for d20 and d09 with the v3 union,
#            using copies of mask_allsites.sh and masked_snps.sh
#            (v2 outputs are left untouched; v3 files get a _v3 / masked_v3 name)
#
# Usage (SNP_call env active, from snp_calling/):
#   bash mask_v3_genets.sh stats
#   bash mask_v3_genets.sh check
#   bash mask_v3_genets.sh depth
#   bash mask_v3_genets.sh paralog
#   bash mask_v3_genets.sh window
#   PARA_WINDOW_BP=<chosen> bash mask_v3_genets.sh paralog   # if the window changes
#   bash mask_v3_genets.sh union
#   bash mask_v3_genets.sh remask          # long: same cost as mask_allsites.sh all
################################################################################

set -euo pipefail

################################################################################
# CONFIGURATION
################################################################################

# Paths: same as paralog_masks.sh / mask_allsites.sh
PROJECT_DIR="/dados04/jorge/CALLICARPA_2"
WORK_DIR="${PROJECT_DIR}/snp_calling"
MASK_DIR="${WORK_DIR}/results/masks"
V3="${MASK_DIR}/v3_genets"                                   # all v3 outputs go here
NATIVE_VCF="${WORK_DIR}/results/native/raw/native13.raw.vcf.gz"  # native-depth joint calls, 13 libraries
D20_SNPS="${WORK_DIR}/results/d20/filtered/cohort.snps.final.vcf.gz" # unmasked d20 SNPs (as in mask_candidates.sh)
HC="${WORK_DIR}/results/het_checks"                          # het_checks.sh output folder
G="${MASK_DIR}/genome.txt"                                   # bedtools genome file (reference order)
GENOME_SIZE=506362408
BT="/home/genomics/miniconda3/envs/paralog_masks/bin/bedtools"
TMP_DIR="${WORK_DIR}/tmp"

# One library per genet (deepest ramet; see header)
GENETS=(547 550 551 552 554 559 561)

# Allele-balance rule (same idea as the "ab" candidate in mask_candidates.sh):
# pooled alt-read fraction across the het calls of the 7 genets outside
# [AB_LO, AB_HI], judged only when the het calls hold >= AB_MIN_READS reads.
# A true heterozygote sits near 0.5; a collapsed paralog where one copy
# carries the alt allele, or copies of unequal number, sits near 0.25 / 0.75.
AB_LO=${AB_LO:-0.3}
AB_HI=${AB_HI:-0.7}
AB_MIN_READS=${AB_MIN_READS:-20}

# All-genets-het rule: a paralog with a fixed difference between two equal
# copies gives balanced reads (~0.5) but makes EVERY tree heterozygous. Flag a
# site where all called genets are het and at least ALLHET_MIN_CALLED are
# called. Set PARA_ALLHET=0 to switch it off.
PARA_ALLHET=${PARA_ALLHET:-1}
ALLHET_MIN_CALLED=${ALLHET_MIN_CALLED:-6}

# Depth layer: a position is high-depth when >= this many of the 7 genets are
# above 2x their own median (same per-tree cutoffs as v2; 7 of 13 ~ half,
# 4 of 7 ~ half).
HIGHDP_MIN_GENETS=${HIGHDP_MIN_GENETS:-4}

# Window around each flagged paralog site (bp each side). Re-chosen in phase
# "window"; rerun phase "paralog" if the chosen value is not the default.
PARA_WINDOW_BP=${PARA_WINDOW_BP:-1000}
WINDOWS_TO_TEST=(0 500 1000)

PHASE="${1:-}"

################################################################################
# HELPERS
################################################################################

log()   { echo "[$(date +'%F %T')] $*"; }
error() { echo "[ERROR] $*" >&2; exit 1; }
bed_bp() { awk '{s+=$3-$2} END{printf "%d", s}' "$1"; }
pct()    { awk -v b="$1" -v g="${GENOME_SIZE}" 'BEGIN{printf "%.2f", 100*b/g}'; }

command -v bcftools >/dev/null || error "conda activate SNP_call"
[[ -x "${BT}" ]] || error "bedtools not found at ${BT}"
mkdir -p "${V3}" "${TMP_DIR}"

# Sample list file, used by bcftools and het_checks.sh
GENET_FILE="${V3}/genets7.txt"
printf '%s\n' "${GENETS[@]}" > "${GENET_FILE}"

STATS="${V3}/site_stats_genets7.tsv.gz"

################################################################################
# PHASE stats — per-site H and allele balance in the 7 genets
################################################################################
# Output columns:
#   chrom pos excesshet13 n_called n_het H het_reads alt_frac D
#   excesshet13 = the original ExcessHet from the 13-library joint calls
#                 (bcftools view -s keeps INFO/ExcessHet as it was), so the
#                 check can compare v2's flag with the genet-level evidence.
#   het_reads   = ref + alt reads summed over the het calls of the 7 genets
#   alt_frac    = alt / het_reads; D = HDplot z (ref - n/2)/sqrt(n/4)
phase_stats() {
    [[ -f "${NATIVE_VCF}.tbi" ]] || error "missing ${NATIVE_VCF}"
    log "stats: per-site H and allele balance in genets ${GENETS[*]}"
    # 1. Keep the 7 genets and biallelic SNPs only
    # 2. Print ExcessHet and GT:AD of each genet
    # 3. awk: count called and het genets; sum ref/alt reads over het calls
    bcftools view -s "$(paste -sd, "${GENET_FILE}")" -m2 -M2 -v snps -Ou "${NATIVE_VCF}" \
      | bcftools query -f '%CHROM\t%POS\t%INFO/ExcessHet[\t%GT:%AD]\n' \
      | awk -F'\t' -v OFS='\t' '
          BEGIN{ print "chrom","pos","excesshet13","n_called","n_het","H","het_reads","alt_frac","D" }
          {
            nc=0; nh=0; r=0; a=0
            for(i=4;i<=NF;i++){
              split($i, x, ":"); g=x[1]
              if(g ~ /\./) continue                      # missing genotype
              nc++
              split(g, al, /[\/|]/)
              if(al[1]!=al[2]){                          # heterozygous call
                nh++
                split(x[2], ad, ",")
                if(ad[1] ~ /^[0-9]+$/) r+=ad[1]
                if(ad[2] ~ /^[0-9]+$/) a+=ad[2]
              }
            }
            n=r+a
            H  = (nc>0) ? sprintf("%.4f", nh/nc) : "NA"
            af = (n>0)  ? sprintf("%.4f", a/n)   : "NA"
            D  = (n>0)  ? sprintf("%.3f", (r-n/2)/sqrt(n/4)) : "NA"
            print $1, $2, $3, nc, nh, H, n, af, D
          }' \
      | bgzip -@ 4 > "${STATS}"
    log "stats: $(zcat "${STATS}" | tail -n +2 | wc -l) biallelic SNPs -> ${STATS}"
}

################################################################################
# PHASE check — what is inside the v2 ExcessHet layer?
################################################################################
# For sites flagged by v2 (ExcessHet >= 20 at 13 libraries) and, as a
# baseline, for all other SNPs, among sites with at least one het genet:
#   evaluable       het calls hold >= AB_MIN_READS reads
#   skewed          evaluable and alt_frac outside [AB_LO, AB_HI]
#   all_het         all called genets het (>= ALLHET_MIN_CALLED called)
#   paralog_like    skewed OR all_het
#   real_like       evaluable, balanced, not all_het -> looks like genuine
#                   heterozygosity, i.e. what v2 should not have masked
# THE NUMBER FOR TATIANA: real_like as % of evaluable, flagged row.
#   high  -> v2 masked a lot of real heterozygosity: rebuild (go on)
#   low   -> the v2 layer is mostly paralogs: the mask can stay
phase_check() {
    [[ -f "${STATS}" ]] || error "run phase stats first"
    local out="${V3}/check_excesshet_layer.tsv"
    log "check: ExcessHet >= 20 sites vs allele balance in 7 genets"
    zcat "${STATS}" | awk -F'\t' -v OFS='\t' \
        -v lo="${AB_LO}" -v hi="${AB_HI}" -v mr="${AB_MIN_READS}" -v ahm="${ALLHET_MIN_CALLED}" '
        NR==1 { next }
        {
          grp = ($3!="." && $3!="" && $3+0>=20) ? "v2_flagged_ExcessHet_ge20" : "v2_not_flagged"
          tot[grp]++
          if($5==0){ nohet[grp]++; next }               # no het genet at this site
          withhet[grp]++
          ev   = ($7>=mr && $8!="NA")
          skew = ev && ($8+0<lo || $8+0>hi)
          ah   = ($4>=ahm && $5==$4)
          if(ev)   nev[grp]++
          if(skew) nsk[grp]++
          if(ah)   nah[grp]++
          if(skew || ah) npl[grp]++
          else if(ev)    nrl[grp]++
        }
        END{
          print "group","snps","no_het_genet","with_het_genet","evaluable","skewed","all_genets_het","paralog_like","real_like","pct_paralog_like_of_evaluable","pct_real_like_of_evaluable"
          split("v2_flagged_ExcessHet_ge20 v2_not_flagged", gs, " ")
          for(k=1;k<=2;k++){ g=gs[k]
            e=nev[g]+0
            print g, tot[g]+0, nohet[g]+0, withhet[g]+0, e, nsk[g]+0, nah[g]+0, npl[g]+0, nrl[g]+0,
                  (e>0? sprintf("%.1f",100*(npl[g]+0)/e) : "NA"),
                  (e>0? sprintf("%.1f",100*(nrl[g]+0)/e) : "NA")
          }
        }' > "${out}"
    printf 'thresholds\tAB_LO=%s AB_HI=%s AB_MIN_READS=%s ALLHET_MIN_CALLED=%s\n' \
        "${AB_LO}" "${AB_HI}" "${AB_MIN_READS}" "${ALLHET_MIN_CALLED}" >> "${out}"
    column -t "${out}"
    log "check: send ${out} to Tatiana before running the next phases"
}

################################################################################
# PHASE depth — high-depth layer counting genets, not libraries
################################################################################
# Reuses the per-tree "above 2x own median" BEDs from paralog_masks.sh
# (results/masks/depth/<tree>.high.bed); no mosdepth rerun.
phase_depth() {
    local out="${V3}/mask_highdepth_genets.bed" beds=()
    for s in "${GENETS[@]}"; do
        [[ -f "${MASK_DIR}/depth/${s}.high.bed" ]] || error "missing ${MASK_DIR}/depth/${s}.high.bed"
        beds+=("${MASK_DIR}/depth/${s}.high.bed")
    done
    log "depth: positions where >= ${HIGHDP_MIN_GENETS} of ${#GENETS[@]} genets are above 2x their median"
    # multiinter column 4 = number of genets whose BED covers the sub-interval
    "${BT}" multiinter -i "${beds[@]}" \
      | awk -v OFS='\t' -v k="${HIGHDP_MIN_GENETS}" '$4>=k {print $1,$2,$3}' \
      | "${BT}" sort -g "${G}" -i - \
      | "${BT}" merge -i - > "${out}"
    local bp; bp=$(bed_bp "${out}")
    log "depth: ${out}: ${bp} bp ($(pct "${bp}") %); v2 layer was $(bed_bp "${MASK_DIR}/mask_highdepth.bed") bp"
}

################################################################################
# PHASE paralog — new paralog sites and windows
################################################################################
paralog_sites() {   # writes ${V3}/mask_paralog_genets_sites.bed
    local out="${V3}/mask_paralog_genets_sites.bed"
    zcat "${STATS}" | awk -F'\t' -v OFS='\t' \
        -v lo="${AB_LO}" -v hi="${AB_HI}" -v mr="${AB_MIN_READS}" \
        -v useah="${PARA_ALLHET}" -v ahm="${ALLHET_MIN_CALLED}" '
        NR==1 || $5==0 { next }                          # header, or no het genet
        {
          skew = ($7>=mr && $8!="NA" && ($8+0<lo || $8+0>hi))
          ah   = (useah==1 && $4>=ahm && $5==$4)
          if(skew || ah) print $1, $2-1, $2
        }' \
      | "${BT}" sort -g "${G}" -i - \
      | "${BT}" merge -i - > "${out}"
    echo "${out}"
}
windows_bed() {     # $1 sites BED, $2 window bp -> stdout
    if [[ "$2" -eq 0 ]]; then cat "$1"
    else "${BT}" slop -b "$2" -g "${G}" -i "$1" | "${BT}" sort -g "${G}" -i - | "${BT}" merge -i -; fi
}
phase_paralog() {
    [[ -f "${STATS}" ]] || error "run phase stats first"
    log "paralog: skewed balance (outside ${AB_LO}-${AB_HI}, >= ${AB_MIN_READS} het reads) OR all genets het (PARA_ALLHET=${PARA_ALLHET})"
    local sites; sites=$(paralog_sites)
    windows_bed "${sites}" "${PARA_WINDOW_BP}" > "${V3}/mask_paralog_genets_windows.bed"
    printf 'AB_LO=%s\nAB_HI=%s\nAB_MIN_READS=%s\nPARA_ALLHET=%s\nALLHET_MIN_CALLED=%s\nPARA_WINDOW_BP=%s\nGENETS=%s\n' \
        "${AB_LO}" "${AB_HI}" "${AB_MIN_READS}" "${PARA_ALLHET}" "${ALLHET_MIN_CALLED}" "${PARA_WINDOW_BP}" "${GENETS[*]}" \
        > "${V3}/thresholds_used.txt"
    local bp; bp=$(bed_bp "${V3}/mask_paralog_genets_windows.bed")
    log "paralog: $(wc -l < "${sites}") site intervals; windows +-${PARA_WINDOW_BP} bp: ${bp} bp ($(pct "${bp}") %); v2 layer was $(bed_bp "${MASK_DIR}/mask_paralog_windows.bed") bp"
}

################################################################################
# PHASE window — choose the window on the 7 genets of d20
################################################################################
# For each window size: union = lowmap + depth by genet + paralog sites
# extended by the window. het_checks.sh is run on unmasked d20 SNPs, 7 genets
# only. Main criterion: mean share of skewed het calls per tree (alt fraction
# < 0.25 or > 0.75 at depth >= 8; true hets give 2-4 %). It does not assume
# independent trees. Fis and the all-het spectrum are reported as secondary:
# the Naguabo genet and its three seedlings are related, so the HWE
# expectation does not strictly apply.
phase_window() {
    local sites="${V3}/mask_paralog_genets_sites.bed" hd="${V3}/mask_highdepth_genets.bed"
    [[ -f "${sites}" ]] || error "run phase paralog first"
    [[ -f "${hd}" ]]    || error "run phase depth first"
    local tab="${V3}/window_choice.tsv"
    printf 'window_bp\tunion_bp\tpct_genome\tsnps_kept\tmean_frac_skewed_het\tall7het_obs_over_exp\tFis_ratio_corrected\texpected_site_Fis\n' > "${tab}"
    for w in "${WINDOWS_TO_TEST[@]}"; do
        local u="${V3}/union_test_w${w}.bed" tag="v3_genets_w${w}"
        cat "${MASK_DIR}/mask_lowmap.bed" "${hd}" <(windows_bed "${sites}" "${w}") \
          | "${BT}" sort -g "${G}" -i - | "${BT}" merge -i - > "${u}"
        log "window ${w} bp: het_checks.sh on d20, 7 genets"
        VCF="${D20_SNPS}" MASK="${u}" TAG="${tag}" SAMPLES="${GENET_FILE}" \
            bash "${WORK_DIR}/het_checks.sh" > "${WORK_DIR}/logs/het_${tag}.log" 2>&1
        local ubp; ubp=$(bed_bp "${u}")
        local f;   f=$(awk -F'\t' -v t="${tag}" '$1==t' "${HC}/fis_summary.tsv")
        local sk;  sk=$(awk -F'\t' 'NR>1 {s+=$6; n++} END{printf "%.4f", s/n}' "${HC}/per_sample_${tag}.tsv")
        local ah;  ah=$(awk -F'\t' -v n="${#GENETS[@]}" '$1==n {print $4}' "${HC}/het_spectrum_${tag}.tsv")
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${w}" "${ubp}" "$(pct "${ubp}")" \
            "$(cut -f3 <<< "${f}")" "${sk}" "${ah}" "$(cut -f8 <<< "${f}")" "$(cut -f9 <<< "${f}")" >> "${tab}"
    done
    column -t "${tab}"
    log "window: choose the smallest window after which mean_frac_skewed_het stops falling; agree it with Tatiana"
}

################################################################################
# PHASE union — v3 union and comparison with v2
################################################################################
phase_union() {
    local hd="${V3}/mask_highdepth_genets.bed" pw="${V3}/mask_paralog_genets_windows.bed"
    local u3="${MASK_DIR}/mask_paralog_union_v3.bed" u2="${MASK_DIR}/mask_paralog_union.bed"
    [[ -f "${hd}" && -f "${pw}" ]] || error "run phases depth and paralog first"
    grep -q "PARA_WINDOW_BP=${PARA_WINDOW_BP}" "${V3}/thresholds_used.txt" \
        || error "paralog windows were built with another PARA_WINDOW_BP; rerun phase paralog"
    cat "${MASK_DIR}/mask_lowmap.bed" "${hd}" "${pw}" \
      | "${BT}" sort -g "${G}" -i - | "${BT}" merge -i - > "${u3}"

    # How much changed: bp only in v2 (now unmasked) and only in v3 (newly masked)
    local s="${V3}/masks_summary_v3.tsv"
    printf 'layer\tv2_bp\tv3_bp\tv2_pct\tv3_pct\tonly_v2_bp\tonly_v3_bp\n' > "${s}"
    cmp_layer() {   # $1 name, $2 v2 bed, $3 v3 bed
        local a b o2 o3
        a=$(bed_bp "$2"); b=$(bed_bp "$3")
        o2=$("${BT}" subtract -a "$2" -b "$3" | awk '{s+=$3-$2} END{printf "%d", s}')
        o3=$("${BT}" subtract -a "$3" -b "$2" | awk '{s+=$3-$2} END{printf "%d", s}')
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "${a}" "${b}" "$(pct "${a}")" "$(pct "${b}")" "${o2}" "${o3}" >> "${s}"
    }
    cmp_layer lowmap    "${MASK_DIR}/mask_lowmap.bed"          "${MASK_DIR}/mask_lowmap.bed"
    cmp_layer highdepth "${MASK_DIR}/mask_highdepth.bed"       "${hd}"
    cmp_layer paralog   "${MASK_DIR}/mask_paralog_windows.bed" "${pw}"
    cmp_layer union     "${u2}"                                "${u3}"
    column -t "${s}"
    log "union: ${u3}"
}

################################################################################
# PHASE remask — masked all-sites and SNP VCFs with the v3 union
################################################################################
# Copies mask_allsites.sh and masked_snps.sh with the v3 mask and v3 output
# names, so the v2 files stay as they are. Same genotype filters (DP < 3,
# het binomial p < 0.01) and max-missing 0.8. Heavy: same run time and disk
# as mask_allsites.sh all (~2.4 GB per dataset).
phase_remask() {
    [[ -f "${MASK_DIR}/mask_paralog_union_v3.bed" ]] || error "run phase union first"
    local a="${TMP_DIR}/mask_allsites_v3.sh" m="${TMP_DIR}/masked_snps_v3.sh"
    sed -e 's/mask_paralog_union\.bed/mask_paralog_union_v3.bed/g' \
        -e 's/masked_v2/masked_v3/g' \
        -e 's/allsites_masked_summary/allsites_masked_v3_summary/g' \
        -e 's/per_tree_heterozygosity/per_tree_heterozygosity_v3/g' \
        -e 's/mask_allsites_${tag}/mask_allsites_v3_${tag}/g' \
        "${WORK_DIR}/mask_allsites.sh" > "${a}"
    sed -e 's/masked_v2/masked_v3/g' "${WORK_DIR}/masked_snps.sh" > "${m}"
    # Safety: refuse to run if a substitution did not take (v2 would be overwritten)
    grep -q 'mask_paralog_union_v3.bed' "${a}"       || error "v3 mask path not set in ${a}"
    grep -q 'cohort.allsites.masked_v3' "${a}"       || error "v3 output name not set in ${a}"
    grep -q 'per_tree_heterozygosity_v3' "${a}"      || error "v3 het table name not set in ${a}"
    grep -q 'cohort.allsites.masked_v3' "${m}"       || error "v3 input not set in ${m}"
    grep -q 'cohort.snps.masked_v3_gtfilt' "${m}"    || error "v3 SNP output not set in ${m}"
    log "remask: all-sites d20 + d09 with mask_paralog_union_v3.bed"
    bash "${a}" all
    for t in d20 d09; do
        log "remask: SNP VCF ${t}"
        bash "${m}" "${t}"
    done
    log "remask: done. Deliver results/<tag>/allsites/cohort.allsites.masked_v3.vcf.gz and results/<tag>/masked/cohort.snps.masked_v3_gtfilt.vcf.gz (d20, d09)"
}

################################################################################
# PHASE reldepth — depth relative to each genet's median (Tatiana, 19 Oct)
################################################################################
# Neither "all genets het" nor allele balance separates paralogs from real
# hets inside the ExcessHet layer (a 2-copy paralog with a fixed difference
# gives AB ~0.5 too). Depth does: a collapsed 2-copy paralog sits near 2x,
# a real het near 1x.
# Per site: FORMAT/DP of each genet / that genet's median (cutoffs.tsv,
# mosdepth on the native BAM); missing DP = 0; site value = mean of the 7.
# GATK FORMAT/DP drops some reads mosdepth counts, so the scale is checked:
# each genet's median FORMAT/DP over not-flagged SNPs is reported, and a
# second normalisation by those VCF medians is computed alongside.
# Groups: not_flagged, flagged (ExcessHet >= 20), flagged_bal_allhet
# (flagged + balanced AB [AB_LO,AB_HI] with >= AB_MIN_READS het reads + all
# called genets het, >= ALLHET_MIN_CALLED called), flagged_rest (flagged minus
# that subset). Also % of each group inside the v2 high-depth layer.
# REGION=Chr17 limits the run (test); outputs then go to reldepth_<REGION>/.
# Outputs (reldepth/): site_reldepth_genets7.tsv.gz (not tracked),
#   reldepth_medians.tsv, reldepth_bins.tsv, reldepth_groups.tsv,
#   reldepth_summary.tsv, reldepth_hist.png/.pdf (plot_reldepth.R)
phase_reldepth() {
    [[ -f "${NATIVE_VCF}.tbi" ]] || error "missing ${NATIVE_VCF}"
    local cut="${MASK_DIR}/depth/cutoffs.tsv" hd="${MASK_DIR}/mask_highdepth.bed"
    local rscript="${RSCRIPT_BIN:-/home/genomics/miniconda3/envs/R_popstat_jorge/bin/Rscript}"
    local region="${REGION:-}" out="${V3}/reldepth"
    [[ -n "${region}" ]] && out="${V3}/reldepth_${region//[:,]/_}"
    mkdir -p "${out}"
    local sites="${out}/site_reldepth_genets7.tsv.gz"

    # Medians in genet order, from cutoffs.tsv
    local meds="" s m
    for s in "${GENETS[@]}"; do
        m=$(awk -v s="${s}" '$1==s{print $2}' "${cut}")
        [[ -n "${m}" ]] || error "no median for ${s} in ${cut}"
        meds+="${m} "
    done
    log "reldepth: medians (cutoffs.tsv) ${GENETS[*]} = ${meds}${region:+ | region ${region}}"

    # 1. Per-site table (BED-like, so it can be intersected with the high-depth layer)
    #    chrom start end flagged bal_allhet dp_<genet>x7 rel_cut
    bcftools view ${region:+-r "${region}"} -s "$(paste -sd, "${GENET_FILE}")" -m2 -M2 -v snps -Ou "${NATIVE_VCF}" \
      | bcftools query -f '%CHROM\t%POS\t%INFO/ExcessHet[\t%GT:%AD:%DP]\n' \
      | awk -F'\t' -v OFS='\t' -v meds="${meds}" -v gen="${GENETS[*]}" \
            -v lo="${AB_LO}" -v hi="${AB_HI}" -v mr="${AB_MIN_READS}" -v ahm="${ALLHET_MIN_CALLED}" '
          BEGIN{ ng=split(meds, md, " "); split(gen, gn, " ")
                 printf "chrom\tstart\tend\tflagged\tbal_allhet"; for(k=1;k<=ng;k++) printf "\tdp_%s", gn[k]; print "\trel_cut" }
          {
            fl = ($3!="." && $3!="" && $3+0>=20)
            nc=0; nh=0; r=0; a=0; rel=0; line=""
            for(i=4;i<=NF;i++){
              split($i, x, ":"); g=x[1]; dp=(x[3] ~ /^[0-9]+$/) ? x[3]+0 : 0
              line = line "\t" dp; rel += dp/md[i-3]
              if(g ~ /\./) continue
              nc++; split(g, al, /[\/|]/)
              if(al[1]!=al[2]){ nh++; split(x[2], ad, ",")
                if(ad[1] ~ /^[0-9]+$/) r+=ad[1]; if(ad[2] ~ /^[0-9]+$/) a+=ad[2] }
            }
            n=r+a
            bal = (n>=mr && a/n>=lo && a/n<=hi)
            ah  = (nc>=ahm && nh==nc)
            printf "%s\t%d\t%d\t%d\t%d%s\t%.4f\n", $1, $2-1, $2, fl, (fl && bal && ah), line, rel/ng
          }' \
      | bgzip -@ 4 > "${sites}"
    log "reldepth: $(zcat "${sites}" | awk 'END{print NR-1}') sites -> ${sites}"

    # 2. Scale check: median FORMAT/DP per genet over not-flagged SNPs
    zcat "${sites}" | awk -F'\t' -v OFS='\t' -v meds="${meds}" -v gen="${GENETS[*]}" '
        BEGIN{ ng=split(meds, md, " "); split(gen, gn, " ") }
        NR>1 && $4==0 { n++; for(k=1;k<=ng;k++){ d=$(5+k); if(d>5000) d=5000; h[k,d]++ } }
        END{ print "genet","median_mosdepth_cutoffs","median_vcf_dp_not_flagged","ratio_vcf_over_mosdepth"
             for(k=1;k<=ng;k++){ c=0; for(d=0; d<=5000; d++){ c+=h[k,d]; if(c>=n/2) break }
               printf "%s\t%s\t%d\t%.3f\n", gn[k], md[k], d, d/md[k] } }' > "${out}/reldepth_medians.tsv"
    column -t "${out}/reldepth_medians.tsv"
    local vmeds
    vmeds=$(awk 'NR>1{printf "%s ", $3}' "${out}/reldepth_medians.tsv")

    # 3. Bins (0.01 wide, 0-6, last bin = >= 6) for both normalisations, plus
    #    per-group counts and high-depth overlap
    zcat "${sites}" | tail -n +2 \
      | "${BT}" intersect -sorted -g "${G}" -c -a - -b "${hd}" \
      | awk -F'\t' -v OFS='\t' -v vmeds="${vmeds}" -v bins="${out}/reldepth_bins.tsv" -v grp="${out}/reldepth_groups.tsv" '
          BEGIN{ ng=split(vmeds, vm, " ") }
          function add(gname, norm, v,   b){ b=int(v/0.01); if(b>600) b=600; H[gname SUBSEP norm SUBSEP b]++ }
          {
            relc=$(6+ng); relv=0
            for(k=1;k<=ng;k++) relv += (vm[k]>0 ? $(5+k)/vm[k] : 0)
            relv/=ng; inhd=($NF>0)
            ng_=1; G1[1]=($4 ? "flagged" : "not_flagged")
            if($4){ ng_=2; G1[2]=($5 ? "flagged_bal_allhet" : "flagged_rest") }
            for(j=1;j<=ng_;j++){ g=G1[j]; N[g]++; HD[g]+=inhd; add(g,"cutoffs",relc); add(g,"vcf",relv) }
          }
          END{
            print "group","norm","bin_lo","n" > bins
            for(key in H){ split(key, p, SUBSEP); printf "%s\t%s\t%.2f\t%d\n", p[1], p[2], p[3]*0.01, H[key] > bins }
            print "group","n_sites","n_in_v2_highdepth" > grp
            for(g in N) print g, N[g], HD[g] > grp
          }'

    # 4. Summary table + histograms
    "${rscript}" "${WORK_DIR}/plot_reldepth.R" "${out}"
    column -t "${out}/reldepth_summary.tsv"
    log "reldepth: done -> ${out}"
}

################################################################################
# DISPATCH (no "all": thresholds and window are agreed between phases)
################################################################################
case "${PHASE}" in
    stats)   phase_stats ;;
    check)   phase_check ;;
    depth)   phase_depth ;;
    paralog) phase_paralog ;;
    window)  phase_window ;;
    union)   phase_union ;;
    remask)  phase_remask ;;
    reldepth) phase_reldepth ;;
    *) error "usage: bash mask_v3_genets.sh stats|check|depth|paralog|window|union|remask|reldepth" ;;
esac
