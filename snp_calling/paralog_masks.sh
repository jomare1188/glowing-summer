#!/bin/bash
################################################################################
# PARALOG MASKS  (feedback.txt step 1)
#
# Both panels (9x and 20x) show genome-wide excess heterozygosity, a deficit of
# homozygous-alt genotypes and very negative per-site Fis. Most likely cause:
# collapsed paralogs, plus C. ampla reads mapped to the C. americana
# reference. --hwe 0.01 cannot remove them with n = 9-13. Three masks, each an
# EXCLUSION BED with the reference's own names (Chr01...):
#
#   mappability  GenMap k=100 e=2 on the reference; mappability < 1 excluded
#   max depth    mosdepth on the 13 native BAMs; a position is excluded where
#                >= HIGHDP_MIN_SAMPLES samples exceed HIGHDP_FACTOR x their
#                own median depth over covered bases
#   paralog      per-SNP ExcessHet, H (het fraction) and D (HDplot allele-
#                balance z) from the native-depth joint calls; thresholds are
#                set after reviewing the diagnostics (phase hetstats)
#
# Why not ExcessHet > 54.69 (GATK default): with n samples, ExcessHet cannot
# exceed the value reached when every sample is heterozygous -- 19.8 (n=9),
# 31.0 (n=13), 45.4 (n=18). A 54.69 cutoff masks nothing here.
#
# Phases:  mappability | depth | hetstats | paralog | union | all
#   hetstats writes stats + diagnostics only; paralog writes the BED using
#   PARA_EXCESSHET_MIN / PARA_H_MIN / PARA_D_MAX (agree on them first).
#
# Environments: paralog_masks (genmap, mosdepth, bedtools) is activated per
# command below; SNP_call must be active for bcftools/tabix.
################################################################################

set -euo pipefail

################################################################################
# CONFIGURATION
################################################################################

PROJECT_DIR="/dados04/jorge/CALLICARPA_2"
REF_GENOME="${PROJECT_DIR}/CallicarpaGenome/car_asm.fa"
WORK_DIR="${PROJECT_DIR}/snp_calling"
GENOME_SIZE=506362408

MASK_DIR="${WORK_DIR}/results/masks"
BAM_DIR="${WORK_DIR}/results/native/bam"
NATIVE_VCF="${WORK_DIR}/results/native/raw/native13.raw.vcf.gz"
MASK_SAMPLES="${WORK_DIR}/results/metrics/samples_d09.txt"
LOG_DIR="${WORK_DIR}/logs"
TMP_DIR="${WORK_DIR}/tmp"
CHECKPOINT_DIR="${WORK_DIR}/checkpoints"
MAIN_LOG="${LOG_DIR}/masks_main.log"

CONDA_BASE="/home/genomics/miniconda3"
MASK_ENV="${CONDA_BASE}/envs/paralog_masks/bin"
RSCRIPT_BIN="${RSCRIPT_BIN:-${CONDA_BASE}/envs/R_popstat_jorge/bin/Rscript}"

THREADS=${THREADS:-16}
DEPTH_JOBS=${DEPTH_JOBS:-7}

# Mappability (Tatiana: GenMap k=100 -e 2)
GENMAP_K=100
GENMAP_E=2

# Max depth: reads counted as in the usable-depth measurement (-F 0xD04 =
# unmapped, secondary, duplicate, supplementary). 0xD04 = 3332.
MOSDEPTH_FLAG=3332
HIGHDP_FACTOR=2
HIGHDP_MIN_SAMPLES=${HIGHDP_MIN_SAMPLES:-7}     # of 13

# Paralog-site thresholds -- placeholders until the diagnostics are reviewed.
PARA_EXCESSHET_MIN=${PARA_EXCESSHET_MIN:-}
PARA_H_MIN=${PARA_H_MIN:-}
PARA_D_MAX=${PARA_D_MAX:-}
PARA_MERGE_BP=${PARA_MERGE_BP:-0}

FORCE_RERUN=${FORCE_RERUN:-false}
PHASE="${1:-all}"

################################################################################
# LOGGING AND CHECKPOINTS  (same helpers as snp_datasets.sh)
################################################################################

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

_emit() { echo -e "$1"; echo -e "$1" | sed 's/\x1b\[[0-9;]*m//g' >> "${MAIN_LOG}"; }
log()   { _emit "${GREEN}[$(date +'%F %T')]${NC} $1"; }
info()  { _emit "${BLUE}[INFO]${NC} $1"; }
warn()  { _emit "${YELLOW}[WARN]${NC} $1"; }
error() { _emit "${RED}[ERROR]${NC} $1"; exit 1; }

checkpoint_exists() { [[ -f "${CHECKPOINT_DIR}/$1.$2.done" ]] && ! $FORCE_RERUN; }
create_checkpoint() { date +%s > "${CHECKPOINT_DIR}/$1.$2.done"; }

# Genome file for bedtools, in .fai order
genome_file() { cut -f1,2 "${REF_GENOME}.fai" > "${MASK_DIR}/genome.txt"; }

bed_bp() { awk '{s+=$3-$2} END{printf "%d", s}' "$1"; }

preflight() {
    mkdir -p "${MASK_DIR}"/{mappability,depth,paralog} "${LOG_DIR}" "${TMP_DIR}" "${CHECKPOINT_DIR}"
    touch "${MAIN_LOG}"
    log "Paralog masks — phase ${PHASE}"
    for t in genmap mosdepth bedtools; do [[ -x "${MASK_ENV}/${t}" ]] || error "${t} missing in ${MASK_ENV}"; done
    for t in bcftools tabix bgzip parallel; do command -v "$t" >/dev/null || error "${t} missing (conda activate SNP_call)"; done
    genome_file
}

################################################################################
# PHASE — MAPPABILITY  (reference only)
################################################################################

phase_mappability() {
    local d="${MASK_DIR}/mappability" idx="${MASK_DIR}/mappability/genmap_index"
    local out="${MASK_DIR}/mask_lowmap.bed"
    if checkpoint_exists masks mappability; then info "  mappability done, skipping"; return 0; fi

    log "  GenMap index"
    rm -rf "${idx}"                          # genmap refuses an existing index dir
    "${MASK_ENV}/genmap" index -F "${REF_GENOME}" -I "${idx}" &> "${LOG_DIR}/masks_genmap_index.log"

    log "  GenMap map -K ${GENMAP_K} -E ${GENMAP_E}"
    "${MASK_ENV}/genmap" map -K "${GENMAP_K}" -E "${GENMAP_E}" -I "${idx}" \
        -O "${d}/car_asm.genmap_k${GENMAP_K}_e${GENMAP_E}" -bg -T "${THREADS}" \
        &> "${LOG_DIR}/masks_genmap_map.log"

    # bedGraph value = 1 / (number of occurrences of the k-mer starting here)
    awk -v OFS='\t' '$4 < 1 {print $1,$2,$3}' "${d}/car_asm.genmap_k${GENMAP_K}_e${GENMAP_E}.bedgraph" \
      | "${MASK_ENV}/bedtools" sort -g "${MASK_DIR}/genome.txt" -i - \
      | "${MASK_ENV}/bedtools" merge -i - > "${out}"
    rm -rf "${idx}"
    log "  mask_lowmap.bed: $(bed_bp "${out}") bp ($(awk -v b="$(bed_bp "${out}")" -v g="${GENOME_SIZE}" 'BEGIN{printf "%.2f", 100*b/g}')% of genome)"
    create_checkpoint masks mappability
}

################################################################################
# PHASE — MAX DEPTH  (native BAMs of the 13 mask samples)
################################################################################

depth_sample() {
    set -euo pipefail
    local s=$1 d="${MASK_DIR}/depth"
    local bam="${BAM_DIR}/${s}.markdup.bam"
    [[ -f "${bam}.bai" ]] || { echo "ERROR: no native BAM for ${s}" >&2; return 1; }
    if checkpoint_exists "masks_${s}" depth; then return 0; fi

    # Pass 1: depth distribution -> median over covered bases
    "${MASK_ENV}/mosdepth" -t 4 -n -F "${MOSDEPTH_FLAG}" "${d}/${s}.pass1" "${bam}"
    # global.dist rows: chrom=total, depth, fraction of bases with >= depth
    local med
    # (array keys are strings in awk: compare as k+0, or "10" sorts before "9")
    med=$(awk '$1=="total" {f[$2+0]=$3} END{c=f[1]; m=0; for(k in f) if(k+0>=1 && f[k]>=c/2 && k+0>m) m=k+0; print m}' \
          "${d}/${s}.pass1.mosdepth.global.dist.txt")
    local cut; cut=$(awk -v m="${med}" -v f="${HIGHDP_FACTOR}" 'BEGIN{printf "%d", m*f}')
    printf '%s\t%s\t%s\n' "${s}" "${med}" "${cut}" > "${d}/${s}.cutoff.tsv"

    # Pass 2: quantize at the cutoff; keep the ">= cut" bin
    MOSDEPTH_Q0=LOW MOSDEPTH_Q1=HIGH \
        "${MASK_ENV}/mosdepth" -t 4 -n -F "${MOSDEPTH_FLAG}" --quantize "0:${cut}:" "${d}/${s}.pass2" "${bam}"
    zcat "${d}/${s}.pass2.quantized.bed.gz" | awk -v OFS='\t' '$4=="HIGH" {print $1,$2,$3}' > "${d}/${s}.high.bed"
    rm -f "${d}/${s}".pass2.quantized.bed.gz* "${d}/${s}".pass1.per-base*
    create_checkpoint "masks_${s}" depth
    log "  depth ${s}: median ${med}x over covered bases, high >= ${cut}x"
}

phase_depth() {
    local out="${MASK_DIR}/mask_highdepth.bed"
    if checkpoint_exists masks highdepth; then info "  highdepth done, skipping"; return 0; fi
    log "  mosdepth on $(wc -l < "${MASK_SAMPLES}") native BAMs"
    PARALLEL_SHELL=/bin/bash parallel --halt soon,fail=1 -j "${DEPTH_JOBS}" depth_sample {} :::: "${MASK_SAMPLES}"

    { printf 'sample\tmedian_covered_X\thigh_cutoff_X\n'; while read -r s; do cat "${MASK_DIR}/depth/${s}.cutoff.tsv"; done < "${MASK_SAMPLES}"; } \
        > "${MASK_DIR}/depth/cutoffs.tsv"

    local beds=(); while read -r s; do beds+=("${MASK_DIR}/depth/${s}.high.bed"); done < "${MASK_SAMPLES}"
    # multiinter column 4 = number of samples covering the sub-interval
    "${MASK_ENV}/bedtools" multiinter -i "${beds[@]}" \
      | awk -v OFS='\t' -v k="${HIGHDP_MIN_SAMPLES}" '$4 >= k {print $1,$2,$3}' \
      | "${MASK_ENV}/bedtools" sort -g "${MASK_DIR}/genome.txt" -i - \
      | "${MASK_ENV}/bedtools" merge -i - > "${out}"
    log "  mask_highdepth.bed: $(bed_bp "${out}") bp ($(awk -v b="$(bed_bp "${out}")" -v g="${GENOME_SIZE}" 'BEGIN{printf "%.2f", 100*b/g}')% of genome), >= ${HIGHDP_MIN_SAMPLES} samples above ${HIGHDP_FACTOR}x median"
    create_checkpoint masks highdepth
}

################################################################################
# PHASE — HET STATS  (diagnostics only; no BED yet)
################################################################################

phase_hetstats() {
    local d="${MASK_DIR}/paralog" stats="${MASK_DIR}/paralog/paralog_site_stats.tsv.gz"
    [[ -f "${NATIVE_VCF}.tbi" ]] || error "native joint VCF missing: ${NATIVE_VCF} (run native_depth.sh)"
    log "  per-site ExcessHet, H, D from $(basename "${NATIVE_VCF}")"
    # Biallelic SNPs only. D = HDplot z-score of the ref-allele read fraction in
    # heterozygotes: (ref - n/2) / sqrt(n/4), with reads summed over hets.
    bcftools view -m2 -M2 -v snps "${NATIVE_VCF}" \
      | bcftools query -f '%CHROM\t%POS\t%INFO/ExcessHet[\t%GT:%AD]\n' \
      | awk -F'\t' -v OFS='\t' 'BEGIN{print "chrom","pos","excesshet","n_called","n_het","H","het_reads","D"}
          { nc=0; nh=0; a=0; b=0
            for(i=4;i<=NF;i++){ split($i,x,":"); g=x[1]; if(g ~ /\./) continue; nc++
              if(g=="0/1"||g=="0|1"||g=="1|0"){ nh++; split(x[2],ad,","); a+=ad[1]; b+=ad[2] } }
            n=a+b; D=(n>0)?(a-n/2)/sqrt(n/4):"NA"; H=(nc>0)?nh/nc:"NA"
            print $1,$2,$3,nc,nh,H,n,D }' \
      | bgzip -@ 4 > "${stats}"

    # Counts at candidate thresholds, for the review with Tatiana
    zcat "${stats}" | awk -F'\t' 'NR>1 {
          n++; if($3>=31) c31++; if($3>=20) c20++; if($3>=13) c13++
          if($6>=0.8) h80++; if($6>=0.6) h60++; if($6>=0.8 && $8!="NA" && ($8>7||$8<-7)) hd++ }
        END{ printf "biallelic_snps\t%d\nexcesshet_ge_31(all_het)\t%d\nexcesshet_ge_20(p<0.01)\t%d\nexcesshet_ge_13(p<0.05)\t%d\nH_ge_0.8\t%d\nH_ge_0.6\t%d\nH_ge_0.8_and_absD_gt_7\t%d\n", n,c31,c20,c13,h80,h60,hd }' \
        > "${d}/candidate_threshold_counts.tsv"
    column -t "${d}/candidate_threshold_counts.tsv" | while read -r l; do info "    ${l}"; done

    "${RSCRIPT_BIN}" - "${stats}" "${d}/paralog_diagnostics" <<'REOF'
args <- commandArgs(TRUE)
x <- read.delim(gzfile(args[1]), na.strings = "NA")
x <- x[is.finite(x$D) & x$n_called >= 10, ]
set.seed(42); s <- x[sample(nrow(x), min(nrow(x), 300000)), ]
draw <- function() {
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
  plot(s$H, s$D, pch = 16, cex = 0.25, col = rgb(0.1, 0.2, 0.5, 0.15),
       xlab = "H (fraction heterozygous)", ylab = "D (allele-balance z)", main = "HDplot, 300K random SNPs")
  abline(h = c(-7, 7), v = c(0.6, 0.8), lty = 2, col = "grey40")
  hist(x$H, breaks = 50, col = "grey70", border = NA, main = "H, all biallelic SNPs", xlab = "H")
  hist(x$excesshet, breaks = 60, col = "grey70", border = NA, main = "ExcessHet (max 31.0 at n = 13)", xlab = "ExcessHet")
  abline(v = c(13, 20), lty = 2, col = "firebrick")
}
png(paste0(args[2], ".png"), width = 1800, height = 600, res = 120); draw(); invisible(dev.off())
pdf(paste0(args[2], ".pdf"), width = 15, height = 5); draw(); invisible(dev.off())
REOF
    log "  diagnostics: ${d}/paralog_diagnostics.png, candidate_threshold_counts.tsv"
    warn "  paralog BED NOT written: agree thresholds, then run phase 'paralog' with PARA_* set"
}

################################################################################
# PHASE — PARALOG BED  (after thresholds are agreed)
################################################################################

phase_paralog() {
    local stats="${MASK_DIR}/paralog/paralog_site_stats.tsv.gz" out="${MASK_DIR}/mask_paralog_sites.bed"
    [[ -f "${stats}" ]] || error "run phase hetstats first"
    [[ -n "${PARA_EXCESSHET_MIN}${PARA_H_MIN}" ]] || error "set PARA_EXCESSHET_MIN and/or PARA_H_MIN (see candidate_threshold_counts.tsv)"
    log "  flag: ExcessHet >= ${PARA_EXCESSHET_MIN:-off} OR (H >= ${PARA_H_MIN:-off} AND |D| > ${PARA_D_MAX:-any}); merge within ${PARA_MERGE_BP} bp"
    zcat "${stats}" | awk -F'\t' -v OFS='\t' -v e="${PARA_EXCESSHET_MIN}" -v h="${PARA_H_MIN}" -v dm="${PARA_D_MAX}" 'NR>1 {
          f=0
          if(e!="" && $3>=e) f=1
          if(h!="" && $6!="NA" && $6>=h && (dm=="" || ($8!="NA" && ($8>dm || $8<-dm)))) f=1
          if(f) print $1,$2-1,$2 }' \
      | "${MASK_ENV}/bedtools" sort -g "${MASK_DIR}/genome.txt" -i - \
      | "${MASK_ENV}/bedtools" merge -d "${PARA_MERGE_BP}" -i - > "${out}"
    log "  mask_paralog_sites.bed: $(wc -l < "${out}") intervals, $(bed_bp "${out}") bp"
    printf 'PARA_EXCESSHET_MIN=%s\nPARA_H_MIN=%s\nPARA_D_MAX=%s\nPARA_MERGE_BP=%s\n' \
        "${PARA_EXCESSHET_MIN}" "${PARA_H_MIN}" "${PARA_D_MAX}" "${PARA_MERGE_BP}" > "${MASK_DIR}/paralog/thresholds_used.txt"
}

################################################################################
# PHASE — UNION, KEEP-SITES AND SUMMARY
################################################################################

# Mean Ho and mean per-site Fis (1 - Ho/He) over biallelic SNPs of a VCF,
# optionally excluding a BED. Prints: n_sites  mean_Ho  mean_Fis
fis_stats() {
    local vcf=$1 excl=${2:-}
    local tflag=(); [[ -n "${excl}" ]] && tflag=(-T "^${excl}")
    bcftools view "${tflag[@]}" -m2 -M2 -v snps "${vcf}" \
      | bcftools query -f '[%GT\t]\n' \
      | awk -F'\t' '{ nc=0; nh=0; alt=0
          for(i=1;i<NF;i++){ g=$i; if(g ~ /\./) continue; nc++; split(g,a,/[\/|]/); alt+=a[1]+a[2]; if(a[1]!=a[2]) nh++ }
          if(nc<2) next; p=alt/(2*nc); he=2*p*(1-p); if(he==0) next
          ho=nh/nc; n++; sho+=ho; sfis+=1-ho/he }
        END{ printf "%d\t%.4f\t%.4f\n", n, sho/n, sfis/n }'
}

phase_union() {
    local union="${MASK_DIR}/mask_paralog_union.bed" keep="${MASK_DIR}/keep_sites.bed"
    local parts=()
    for m in mask_lowmap mask_highdepth mask_paralog_sites; do
        [[ -f "${MASK_DIR}/${m}.bed" ]] && parts+=("${MASK_DIR}/${m}.bed") || warn "  ${m}.bed missing -- union built without it"
    done
    (( ${#parts[@]} > 0 )) || error "no masks built yet"
    cat "${parts[@]}" | "${MASK_ENV}/bedtools" sort -g "${MASK_DIR}/genome.txt" -i - \
      | "${MASK_ENV}/bedtools" merge -i - > "${union}"
    "${MASK_ENV}/bedtools" complement -i "${union}" -g "${MASK_DIR}/genome.txt" > "${keep}"

    local summ="${MASK_DIR}/masks_summary.tsv" R="${WORK_DIR}/results"
    printf 'mask\tbp\tpct_genome\td09_final_snps_removed\td20_final_snps_removed\n' > "${summ}"
    for b in "${parts[@]}" "${union}"; do
        local bp; bp=$(bed_bp "${b}")
        local r09 r20
        r09=$(bcftools view -H -T "${b}" "${R}/d09/filtered/cohort.snps.final.vcf.gz" | wc -l)
        r20=$(bcftools view -H -T "${b}" "${R}/d20/filtered/cohort.snps.final.vcf.gz" | wc -l)
        printf '%s\t%s\t%.2f\t%s\t%s\n' "$(basename "${b}" .bed)" "${bp}" \
            "$(awk -v b="${bp}" -v g="${GENOME_SIZE}" 'BEGIN{print 100*b/g}')" "${r09}" "${r20}" >> "${summ}"
    done
    column -t "${summ}" | while read -r l; do info "    ${l}"; done

    # Does masking fix the excess heterozygosity? Mean Ho / Fis before vs after.
    local fis="${MASK_DIR}/fis_before_after.tsv"
    printf 'dataset\tset\tn_snps\tmean_Ho\tmean_Fis\n' > "${fis}"
    for t in d09 d20; do
        local v="${R}/${t}/filtered/cohort.snps.final.vcf.gz"
        printf '%s\tunmasked\t%s\n' "${t}" "$(fis_stats "${v}")" >> "${fis}"
        printf '%s\tmasked\t%s\n'   "${t}" "$(fis_stats "${v}" "${union}")" >> "${fis}"
    done
    column -t "${fis}" | while read -r l; do info "    ${l}"; done
    log "  union: ${union}; keep: ${keep}"
}

export -f depth_sample _emit log info warn error checkpoint_exists create_checkpoint
export MASK_DIR BAM_DIR MASK_ENV MOSDEPTH_FLAG HIGHDP_FACTOR LOG_DIR MAIN_LOG CHECKPOINT_DIR \
       FORCE_RERUN RED GREEN YELLOW BLUE NC

main() {
    preflight
    case "${PHASE}" in
        mappability) phase_mappability ;;
        depth)       phase_depth ;;
        hetstats)    phase_hetstats ;;
        paralog)     phase_paralog ;;
        union)       phase_union ;;
        all)         phase_mappability; phase_depth; phase_hetstats ;;
        *)           error "unknown phase '${PHASE}' (mappability|depth|hetstats|paralog|union|all)" ;;
    esac
    log "paralog_masks.sh finished (phase ${PHASE})"
}

main
