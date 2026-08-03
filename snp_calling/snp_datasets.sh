#!/bin/bash
################################################################################
# SNP CALLING — DEPTH-EQUALISED DATASETS (9x and 20x)
#
# Rationale: mapping rates in this cohort range from 6% to 77% against nearly
# uniform raw sequencing depth, so realised coverage varies ~19-fold between
# samples. Depth variance of that magnitude biases which SNPs are callable and
# how genotypes are assigned. This pipeline removes the variance by building two
# datasets in which EVERY sample is downsampled to the same depth.
#
#   9x  dataset: samples reaching >= 9x,  each downsampled to exactly  9x
#   20x dataset: samples reaching >= 20x, each downsampled to exactly 20x
#
# Per dataset it produces: a filtered SNP set, an LD-pruned set, 50K and 5K
# random subsets, and the MAF 0.05 / 0.10 inputs for ROH analysis.
#
# Storage discipline (kept from snp_optimized.sh and extended):
#   - alignment is streamed straight into sorting, no unsorted BAM ever lands
#   - the sorted BAM is deleted as soon as duplicates are marked
#   - both downsampled BAMs are cut from ONE markdup BAM, then deleted after
#     their GVCF exists
#   - the markdup BAM is deleted once every dataset has consumed it
#   - every step is checkpointed, so an interrupted run resumes where it stopped
#
# Usage:  bash snp_datasets.sh [phase]
#           phase = all | align | assign | call | downstream   (default: all)
################################################################################

set -euo pipefail

################################################################################
# CONFIGURATION
################################################################################

PROJECT_DIR="/dados04/jorge/CALLICARPA_2"
READS_DIR="${PROJECT_DIR}/qc/qc_results/trimmed_reads"
REF_GENOME="${PROJECT_DIR}/CallicarpaGenome/car_asm.fa"
WORK_DIR="${PROJECT_DIR}/snp_calling"

OUT_DIR="${WORK_DIR}/results"
LOG_DIR="${WORK_DIR}/logs"
TMP_DIR="${WORK_DIR}/tmp"
CHECKPOINT_DIR="${WORK_DIR}/checkpoints"
MAIN_LOG="${LOG_DIR}/pipeline_main.log"

GENOME_SIZE=506362408

# Datasets to build: <tag>:<target depth in X>
DATASETS=("d09:9" "d20:20")

# Resources
TOTAL_CORES=${TOTAL_CORES:-64}
MAX_PARALLEL_SAMPLES=${MAX_PARALLEL_SAMPLES:-6}
THREADS_PER_SAMPLE=${THREADS_PER_SAMPLE:-10}
MEM_PER_SAMPLE="${MEM_PER_SAMPLE:-40G}"
MEM_COMBINE="${MEM_COMBINE:-100G}"
MEM_GENOTYPE="${MEM_GENOTYPE:-200G}"

# Quality thresholds
MIN_BQ=20
SUBSAMPLE_SEED=42          # fixed so downsampling is reproducible

# GATK hard filters (site level, depth independent)
HARD_FILTERS=(
    "QD_filter:QD < 2.0"
    "FS_filter:FS > 60.0"
    "MQ_filter:MQ < 40.0"
    "SOR_filter:SOR > 3.0"
    "MQRankSum_filter:MQRankSum < -12.5"
    "ReadPosRankSum_filter:ReadPosRankSum < -8.0"
)

# Population filters. Depth bounds are NOT constants here: they are derived from
# each dataset's target depth further down. See scale_depth_filters().
POP_MAF=0.01
POP_HWE=0.01
POP_MAX_MISSING=0.8

# ROH input filters (README: no HWE filter, it removes the signal being measured)
ROH_MAX_MISSING=0.9
ROH_MAFS=(0.05 0.10)

# LD pruning
LD_WINDOW=50
LD_STEP=10
LD_R2=0.2
SUBSET_SIZES=(50000 5000)

# Storage
KEEP_SORTED_BAM=false      # delete after MarkDuplicates
KEEP_MARKDUP_BAM=false     # delete once every dataset has been downsampled from it
KEEP_DOWNSAMPLED_BAM=false # delete once the GVCF exists
KEEP_COMBINED_GVCF=false   # delete after GenotypeGVCFs
DISK_WARN_PCT=85
DISK_ABORT_PCT=96

FORCE_RERUN=${FORCE_RERUN:-false}

PHASE="${1:-all}"

################################################################################
# LOGGING AND CHECKPOINTS
################################################################################

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; BLUE='\033[0;34m'; NC='\033[0m'

_emit() { echo -e "$1"; echo -e "$1" | sed 's/\x1b\[[0-9;]*m//g' >> "${MAIN_LOG}"; }
log()   { _emit "${GREEN}[$(date +'%F %T')]${NC} $1"; }
info()  { _emit "${BLUE}[INFO]${NC} $1"; }
warn()  { _emit "${YELLOW}[WARN]${NC} $1"; }
error() { _emit "${RED}[ERROR]${NC} $1"; exit 1; }

checkpoint_exists() { [[ -f "${CHECKPOINT_DIR}/$1.$2.done" ]] && ! $FORCE_RERUN; }
create_checkpoint() { date +%s > "${CHECKPOINT_DIR}/$1.$2.done"; }

disk_pct() { df -P "${WORK_DIR}" | awk 'NR==2 {gsub(/%/,"",$5); print $5}'; }

check_disk() {
    local pct; pct=$(disk_pct)
    (( pct >= DISK_ABORT_PCT )) && error "disk ${pct}% full, aborting to avoid a truncated write"
    (( pct >= DISK_WARN_PCT )) && warn "disk ${pct}% full"
    return 0
}

################################################################################
# PHASE 0 — PREFLIGHT
################################################################################

preflight() {
    mkdir -p "${OUT_DIR}"/{alignment,markdup,metrics,downsampled} \
             "${LOG_DIR}" "${TMP_DIR}" "${CHECKPOINT_DIR}"
    : > "${MAIN_LOG}" 2>/dev/null || true

    log "Preflight"

    local missing=()
    for t in bwa samtools gatk bcftools tabix bgzip parallel plink plink2 vcftools; do
        command -v "$t" >/dev/null || missing+=("$t")
    done
    if (( ${#missing[@]} > 0 )); then
        error "missing tools: ${missing[*]}
  Activate both environments:
    source /home/genomics/miniconda3/etc/profile.d/conda.sh
    conda activate SNP_call
    export PATH=\"/home/genomics/miniconda3/envs/popgen_tools/bin:\$PATH\""
    fi

    [[ -f "${REF_GENOME}" ]] || error "reference not found: ${REF_GENOME}"
    [[ -d "${READS_DIR}" ]]  || error "trimmed reads not found: ${READS_DIR} (run qc/qc_pipeline_v2.sh first)"

    # Reference indexes
    [[ -f "${REF_GENOME}.bwt" ]] || { log "  building BWA index"; bwa index "${REF_GENOME}" &> "${LOG_DIR}/bwa_index.log"; }
    [[ -f "${REF_GENOME}.fai" ]] || { log "  building FASTA index"; samtools faidx "${REF_GENOME}"; }
    local dict="${REF_GENOME%.fa}.dict"
    [[ -f "${dict}" ]] || { log "  building sequence dictionary"; gatk CreateSequenceDictionary -R "${REF_GENOME}" -O "${dict}" &> "${LOG_DIR}/dict.log"; }

    # Sample discovery, both mates required
    SAMPLES=()
    shopt -s nullglob
    for r1 in "${READS_DIR}"/*_R1_trimmed.fastq.gz; do
        local s r2
        s=$(basename "${r1}" _R1_trimmed.fastq.gz)
        r2="${READS_DIR}/${s}_R2_trimmed.fastq.gz"
        [[ -f "${r2}" ]] || error "missing mate for ${r1}"
        SAMPLES+=("${s}")
    done
    shopt -u nullglob
    (( ${#SAMPLES[@]} > 0 )) || error "no *_R1_trimmed.fastq.gz in ${READS_DIR}"

    log "  ${#SAMPLES[@]} samples: ${SAMPLES[*]}"
    log "  disk $(disk_pct)% used, $(df -Ph "${WORK_DIR}" | awk 'NR==2{print $4}') free"
    check_disk
}

################################################################################
# PHASE 1 — ALIGN, MARK DUPLICATES, MEASURE DEPTH
################################################################################
# Runs once at full depth. Both datasets are cut from these BAMs, so alignment
# is never repeated per target depth.

process_sample() {
    set -euo pipefail
    local sample=$1
    local r1="${READS_DIR}/${sample}_R1_trimmed.fastq.gz"
    local r2="${READS_DIR}/${sample}_R2_trimmed.fastq.gz"
    local sorted="${OUT_DIR}/alignment/${sample}.sorted.bam"
    local markdup="${OUT_DIR}/markdup/${sample}.markdup.bam"
    local depth_file="${OUT_DIR}/metrics/${sample}.depth.txt"

    # -- alignment ----------------------------------------------------------
    if checkpoint_exists "${sample}" alignment; then
        info "  ${sample}: alignment done, skipping"
    else
        log "  ${sample}: aligning"
        check_disk
        # Streamed into sort: an unsorted BAM is never written to disk.
        bwa mem -t "${THREADS_PER_SAMPLE}" -M \
            -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}\tPU:unit1" \
            "${REF_GENOME}" "${r1}" "${r2}" 2> "${LOG_DIR}/${sample}_bwa.log" \
          | samtools sort -@ 4 -m 2G -l 9 -T "${TMP_DIR}/${sample}.sort" -o "${sorted}" - \
            2> "${LOG_DIR}/${sample}_sort.log"
        samtools index -@ 4 "${sorted}"
        samtools flagstat -@ 4 "${sorted}" > "${OUT_DIR}/metrics/${sample}_flagstat.txt"
        create_checkpoint "${sample}" alignment
    fi

    # -- duplicate marking --------------------------------------------------
    if checkpoint_exists "${sample}" markdup; then
        info "  ${sample}: markdup done, skipping"
    else
        log "  ${sample}: marking duplicates"
        check_disk
        mkdir -p "${TMP_DIR}/${sample}"
        gatk --java-options "-Xmx${MEM_PER_SAMPLE}" MarkDuplicates \
            -I "${sorted}" -O "${markdup}" \
            -M "${OUT_DIR}/metrics/${sample}_duplicate_metrics.txt" \
            --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT \
            --COMPRESSION_LEVEL 9 --TMP_DIR "${TMP_DIR}/${sample}" \
            &> "${LOG_DIR}/${sample}_markdup.log"
        rm -rf "${TMP_DIR}/${sample}"
        create_checkpoint "${sample}" markdup
        if ! ${KEEP_SORTED_BAM}; then
            rm -f "${sorted}" "${sorted}.bai"
            info "  ${sample}: sorted BAM removed"
        fi
    fi

    # -- depth measurement --------------------------------------------------
    # "bases mapped (cigar)" counts aligned bases only, so soft-clipped tails do
    # not inflate the estimate. Duplicates, secondary and supplementary
    # alignments are excluded (-F 0xD04) because none of them contribute
    # independent evidence to genotyping.
    if checkpoint_exists "${sample}" depth; then
        info "  ${sample}: depth done, skipping"
    else
        log "  ${sample}: measuring depth"
        local bases
        # No early `exit` in awk: closing the pipe would SIGPIPE samtools, and
        # pipefail would then fail the whole step on an otherwise good BAM.
        bases=$(samtools stats -@ 4 -F 0xD04 "${markdup}" \
                | awk -F'\t' '/^SN\tbases mapped \(cigar\)/ {v=$3} END{print v}')
        [[ -n "${bases}" ]] || { echo "ERROR: ${sample}: could not read mapped bases" >&2; return 1; }
        awk -v b="${bases}" -v g="${GENOME_SIZE}" 'BEGIN{printf "%.4f\n", b/g}' > "${depth_file}"
        log "  ${sample}: $(cat "${depth_file}")x usable depth"
        create_checkpoint "${sample}" depth
    fi
}
export -f process_sample

phase_align() {
    log "PHASE 1 — alignment, duplicate marking, depth measurement"
    export READS_DIR OUT_DIR LOG_DIR TMP_DIR CHECKPOINT_DIR REF_GENOME GENOME_SIZE
    export THREADS_PER_SAMPLE MEM_PER_SAMPLE KEEP_SORTED_BAM FORCE_RERUN MAIN_LOG
    export RED GREEN YELLOW BLUE NC WORK_DIR DISK_WARN_PCT DISK_ABORT_PCT
    export -f log info warn error _emit checkpoint_exists create_checkpoint disk_pct check_disk

    printf '%s\n' "${SAMPLES[@]}" \
      | parallel --halt soon,fail=1 -j "${MAX_PARALLEL_SAMPLES}" process_sample {}
    log "PHASE 1 complete"
}

################################################################################
# PHASE 2 — DATASET ASSIGNMENT
################################################################################
# A sample joins a dataset only if its measured depth reaches that target.
# Downsampling can remove reads, never invent them.

phase_assign() {
    log "PHASE 2 — assigning samples to datasets"
    local table="${OUT_DIR}/metrics/depth_summary.tsv"
    printf 'sample\tusable_depth_X' > "${table}"
    for d in "${DATASETS[@]}"; do printf '\t%s' "${d%%:*}" >> "${table}"; done
    printf '\n' >> "${table}"

    for tag_target in "${DATASETS[@]}"; do
        : > "${OUT_DIR}/metrics/samples_${tag_target%%:*}.txt"
    done

    for s in "${SAMPLES[@]}"; do
        local df="${OUT_DIR}/metrics/${s}.depth.txt"
        [[ -f "${df}" ]] || error "depth missing for ${s}; run the align phase first"
        local depth; depth=$(cat "${df}")
        printf '%s\t%s' "${s}" "${depth}" >> "${table}"
        for tag_target in "${DATASETS[@]}"; do
            local tag="${tag_target%%:*}" target="${tag_target##*:}"
            if awk -v d="${depth}" -v t="${target}" 'BEGIN{exit !(d>=t)}'; then
                echo "${s}" >> "${OUT_DIR}/metrics/samples_${tag}.txt"
                printf '\tIN' >> "${table}"
            else
                printf '\tout' >> "${table}"
            fi
        done
        printf '\n' >> "${table}"
    done

    log "  depth summary written to ${table}"
    column -t "${table}" | sed 's/^/    /' | while read -r l; do info "${l}"; done
    for tag_target in "${DATASETS[@]}"; do
        local tag="${tag_target%%:*}" target="${tag_target##*:}"
        local n; n=$(wc -l < "${OUT_DIR}/metrics/samples_${tag}.txt")
        log "  dataset ${tag} (${target}x): ${n} samples"
        (( n >= 2 )) || error "dataset ${tag} has ${n} samples, too few to genotype"
    done
    create_checkpoint cohort assignment
}

################################################################################
# PHASE 3a — DOWNSAMPLE AND CALL PER DATASET
################################################################################

downsample_and_call() {
    set -euo pipefail
    local tag=$1 target=$2 sample=$3
    local markdup="${OUT_DIR}/markdup/${sample}.markdup.bam"
    local ds_bam="${OUT_DIR}/downsampled/${tag}/${sample}.${tag}.bam"
    local gvcf="${OUT_DIR}/${tag}/gvcf/${sample}.g.vcf.gz"

    if checkpoint_exists "${sample}" "${tag}_gvcf"; then
        info "  ${tag}/${sample}: GVCF done, skipping"
        return 0
    fi

    local depth frac
    depth=$(cat "${OUT_DIR}/metrics/${sample}.depth.txt")
    frac=$(awk -v t="${target}" -v d="${depth}" 'BEGIN{f=t/d; if(f>1)f=1; printf "%.5f", f}')

    if ! checkpoint_exists "${sample}" "${tag}_downsample"; then
        log "  ${tag}/${sample}: ${depth}x -> ${target}x (keeping ${frac} of reads)"
        check_disk
        mkdir -p "$(dirname "${ds_bam}")"
        if awk -v f="${frac}" 'BEGIN{exit !(f>=0.999)}'; then
            # Already at target; hard-link rather than copy a second BAM to disk.
            ln -f "${markdup}" "${ds_bam}"
            ln -f "${markdup}.bai" "${ds_bam}.bai" 2>/dev/null || samtools index -@ 4 "${ds_bam}"
        else
            # --subsample hashes read NAMES, so both mates are kept or dropped
            # together and the pairing stays intact.
            samtools view -@ 4 -b --subsample "${frac}" --subsample-seed "${SUBSAMPLE_SEED}" \
                -o "${ds_bam}" "${markdup}"
            samtools index -@ 4 "${ds_bam}"
        fi
        create_checkpoint "${sample}" "${tag}_downsample"
    fi

    log "  ${tag}/${sample}: HaplotypeCaller"
    check_disk
    mkdir -p "$(dirname "${gvcf}")" "${TMP_DIR}/${tag}_${sample}"
    gatk --java-options "-Xmx${MEM_PER_SAMPLE}" HaplotypeCaller \
        -R "${REF_GENOME}" -I "${ds_bam}" -O "${gvcf}" \
        -ERC GVCF \
        --native-pair-hmm-threads "${THREADS_PER_SAMPLE}" \
        --min-base-quality-score "${MIN_BQ}" \
        --tmp-dir "${TMP_DIR}/${tag}_${sample}" \
        &> "${LOG_DIR}/${sample}_${tag}_hc.log"
    rm -rf "${TMP_DIR}/${tag}_${sample}"
    create_checkpoint "${sample}" "${tag}_gvcf"

    if ! ${KEEP_DOWNSAMPLED_BAM}; then
        rm -f "${ds_bam}" "${ds_bam}.bai"
        info "  ${tag}/${sample}: downsampled BAM removed"
    fi
}
export -f downsample_and_call

################################################################################
# PHASE 3b — JOINT GENOTYPING AND FILTERING PER DATASET
################################################################################

# Depth-based population filters must track the dataset's depth. Copying the
# original --min-meanDP 20 onto the 9x dataset would discard every site: mean
# depth there is 9 by construction. Bounds are set to half and three times the
# target, which keeps well-behaved sites while still removing collapsed repeats.
scale_depth_filters() {
    local target=$1
    MIN_MEAN_DP=$(awk -v t="${target}" 'BEGIN{printf "%d", (t/2 < 3 ? 3 : t/2)}')
    MAX_MEAN_DP=$(awk -v t="${target}" 'BEGIN{printf "%d", t*3}')
}

genotype_dataset() {
    local tag=$1 target=$2
    local dd="${OUT_DIR}/${tag}"
    mkdir -p "${dd}"/{gvcf,raw,filtered,downstream,roh}

    local combined="${dd}/raw/cohort.g.vcf.gz"
    local raw="${dd}/raw/cohort.raw.vcf.gz"

    if ! checkpoint_exists "${tag}" joint_genotyping; then
        log "  ${tag}: combining GVCFs"
        local args=()
        while read -r s; do args+=(-V "${dd}/gvcf/${s}.g.vcf.gz"); done < "${OUT_DIR}/metrics/samples_${tag}.txt"
        [[ -f "${combined}" ]] || gatk --java-options "-Xmx${MEM_COMBINE}" CombineGVCFs \
            -R "${REF_GENOME}" "${args[@]}" -O "${combined}" &> "${LOG_DIR}/${tag}_combine.log"

        log "  ${tag}: joint genotyping"
        [[ -f "${raw}" ]] || gatk --java-options "-Xmx${MEM_GENOTYPE}" GenotypeGVCFs \
            -R "${REF_GENOME}" -V "${combined}" -O "${raw}" &> "${LOG_DIR}/${tag}_genotype.log"
        create_checkpoint "${tag}" joint_genotyping

        if ! ${KEEP_COMBINED_GVCF}; then
            rm -f "${combined}" "${combined}.tbi"
            info "  ${tag}: combined GVCF removed (per-sample GVCFs retained)"
        fi
    else
        info "  ${tag}: joint genotyping done, skipping"
    fi

    if checkpoint_exists "${tag}" filtering; then
        info "  ${tag}: filtering done, skipping"
        return 0
    fi

    scale_depth_filters "${target}"
    log "  ${tag}: depth filters scaled to target ${target}x -> meanDP ${MIN_MEAN_DP}-${MAX_MEAN_DP}"

    local snps="${dd}/filtered/cohort.snps.vcf.gz"
    local marked="${dd}/filtered/cohort.snps.marked.vcf.gz"
    local pass="${dd}/filtered/cohort.snps.pass.vcf.gz"
    local final="${dd}/filtered/cohort.snps.final.vcf.gz"

    log "  ${tag}: selecting SNPs"
    [[ -f "${snps}" ]] || gatk SelectVariants -R "${REF_GENOME}" -V "${raw}" \
        -select-type SNP -O "${snps}" &> "${LOG_DIR}/${tag}_select.log"

    log "  ${tag}: GATK hard filters"
    if [[ ! -f "${marked}" ]]; then
        local fargs=()
        for f in "${HARD_FILTERS[@]}"; do fargs+=(--filter-name "${f%%:*}" --filter-expression "${f#*:}"); done
        gatk VariantFiltration -R "${REF_GENOME}" -V "${snps}" -O "${marked}" \
            "${fargs[@]}" &> "${LOG_DIR}/${tag}_filtration.log"
    fi

    log "  ${tag}: extracting PASS"
    if [[ ! -f "${pass}" ]]; then
        bcftools view -f PASS -O z -o "${pass}" "${marked}"
        tabix -f -p vcf "${pass}"
    fi

    log "  ${tag}: population filters (maf ${POP_MAF}, hwe ${POP_HWE}, max-missing ${POP_MAX_MISSING}, meanDP ${MIN_MEAN_DP}-${MAX_MEAN_DP})"
    if [[ ! -f "${final}" ]]; then
        vcftools --gzvcf "${pass}" \
            --remove-indels \
            --hwe "${POP_HWE}" \
            --maf "${POP_MAF}" \
            --max-missing "${POP_MAX_MISSING}" \
            --min-meanDP "${MIN_MEAN_DP}" \
            --max-meanDP "${MAX_MEAN_DP}" \
            --recode --recode-INFO-all --stdout 2> "${LOG_DIR}/${tag}_popfilter.log" \
          | bgzip -@ 4 > "${final}"
        tabix -f -p vcf "${final}"
    fi

    bcftools stats "${final}" > "${dd}/filtered/cohort.snps.final.stats.txt"
    log "  ${tag}: $(bcftools index -n "${final}") SNPs in final set"

    # Intermediates are regenerable from cohort.raw.vcf.gz; drop them.
    rm -f "${snps}" "${snps}.tbi" "${marked}" "${marked}.tbi"
    create_checkpoint "${tag}" filtering
}

################################################################################
# PHASE 3c — DOWNSTREAM SETS PER DATASET
################################################################################

ld_prune_and_subset() {
    local tag=$1
    local dd="${OUT_DIR}/${tag}"
    local final="${dd}/filtered/cohort.snps.final.vcf.gz"
    local ddir="${dd}/downstream"
    local pruned="${ddir}/ld_pruned.snps.vcf.gz"

    if checkpoint_exists "${tag}" ld_pruning; then
        info "  ${tag}: LD pruning done, skipping"
    else
        log "  ${tag}: LD pruning (window ${LD_WINDOW}, step ${LD_STEP}, r2 ${LD_R2})"
        plink --vcf "${final}" --make-bed --out "${ddir}/tmp_data" \
              --allow-extra-chr --double-id --set-missing-var-ids '@:#' \
              &> "${LOG_DIR}/${tag}_plink_import.log"
        plink --bfile "${ddir}/tmp_data" \
              --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
              --allow-extra-chr --out "${ddir}/ld_pruning" \
              &> "${LOG_DIR}/${tag}_plink_prune.log"
        plink --bfile "${ddir}/tmp_data" --extract "${ddir}/ld_pruning.prune.in" \
              --recode vcf-iid bgz --allow-extra-chr --out "${ddir}/ld_pruned.snps" \
              &> "${LOG_DIR}/${tag}_plink_extract.log"
        tabix -f -p vcf "${pruned}"
        bcftools stats "${pruned}" > "${ddir}/ld_pruned.snps.stats.txt"
        rm -f "${ddir}"/tmp_data.* "${ddir}"/ld_pruning.{prune.in,prune.out,log,nosex}
        log "  ${tag}: $(bcftools index -n "${pruned}") SNPs after LD pruning"
        create_checkpoint "${tag}" ld_pruning
    fi

    # Random subsets of the pruned panel
    local total; total=$(bcftools index -n "${pruned}")
    for n in "${SUBSET_SIZES[@]}"; do
        local label out
        label=$(( n / 1000 ))K
        out="${ddir}/ld_pruned_${label}.snps.vcf.gz"
        if checkpoint_exists "${tag}" "subset_${label}"; then
            info "  ${tag}: ${label} subset done, skipping"
            continue
        fi
        if (( total < n )); then
            warn "  ${tag}: only ${total} pruned SNPs available, cannot draw ${n}; writing full pruned set as ${label}"
            cp -f "${pruned}" "${out}"; cp -f "${pruned}.tbi" "${out}.tbi"
        else
            log "  ${tag}: drawing ${n} random SNPs"
            # shuf needs the body only; the header is prepended untouched so the
            # result stays a valid VCF, and re-sorting keeps it tabix-indexable.
            { bcftools view -h "${pruned}"
              bcftools view -H "${pruned}" | shuf -n "${n}" --random-source=<(yes "${SUBSAMPLE_SEED}") \
                | sort -k1,1 -k2,2n
            } | bgzip -@ 4 > "${out}"
            tabix -f -p vcf "${out}"
        fi
        bcftools stats "${out}" > "${out%.vcf.gz}.stats.txt"
        create_checkpoint "${tag}" "subset_${label}"
    done
}

prepare_roh_input() {
    local tag=$1
    local dd="${OUT_DIR}/${tag}"
    # Deliberately the FULL filtered SNP set, not the LD-pruned panel: pruning
    # removes the correlated markers that ROH detection needs to trace
    # contiguous homozygous tracts, and the 5K/50K subsets are far too sparse
    # for segment boundaries on a 506 Mb genome.
    local final="${dd}/filtered/cohort.snps.final.vcf.gz"
    local rdir="${dd}/roh"

    if checkpoint_exists "${tag}" roh_input; then
        info "  ${tag}: ROH inputs done, skipping"
        return 0
    fi

    local biallelic="${rdir}/biallelic.tmp.vcf.gz"
    log "  ${tag}: extracting biallelic SNPs for ROH"
    bcftools view --min-alleles 2 --max-alleles 2 --type snps \
        -O z -o "${biallelic}" "${final}"
    tabix -f -p vcf "${biallelic}"

    for maf in "${ROH_MAFS[@]}"; do
        local mtag="maf${maf/./}"
        local out="${rdir}/for_roh_dense_${mtag}.vcf.gz"
        # No HWE filter here: inbred individuals violate HWE at exactly the
        # heterozygous loci that carry the ROH signal.
        log "  ${tag}: ROH input MAF ${maf} (no HWE filter)"
        vcftools --gzvcf "${biallelic}" \
            --max-missing "${ROH_MAX_MISSING}" --maf "${maf}" \
            --recode --recode-INFO-all --stdout 2> "${LOG_DIR}/${tag}_roh_${mtag}.log" \
          | bgzip -@ 4 > "${out}"
        tabix -f -p vcf "${out}"

        plink2 --vcf "${out}" --make-bed --out "${rdir}/callicarpa_${mtag}" \
               --allow-extra-chr --max-alleles 2 \
               --set-missing-var-ids '@:#' --new-id-max-allele-len 50 \
               &> "${LOG_DIR}/${tag}_plink2_${mtag}.log"
        log "  ${tag}: ${mtag} -> $(bcftools index -n "${out}") SNPs"
    done
    rm -f "${biallelic}" "${biallelic}.tbi"
    create_checkpoint "${tag}" roh_input
}

################################################################################
# PHASE ORCHESTRATION
################################################################################

phase_call() {
    log "PHASE 3 — downsampling, calling and joint genotyping"
    export OUT_DIR LOG_DIR TMP_DIR CHECKPOINT_DIR REF_GENOME MAIN_LOG WORK_DIR
    export THREADS_PER_SAMPLE MEM_PER_SAMPLE MIN_BQ SUBSAMPLE_SEED
    export KEEP_DOWNSAMPLED_BAM FORCE_RERUN DISK_WARN_PCT DISK_ABORT_PCT
    export RED GREEN YELLOW BLUE NC
    export -f log info warn error _emit checkpoint_exists create_checkpoint disk_pct check_disk

    for tag_target in "${DATASETS[@]}"; do
        local tag="${tag_target%%:*}" target="${tag_target##*:}"
        log "  dataset ${tag} (target ${target}x)"
        parallel --halt soon,fail=1 -j "${MAX_PARALLEL_SAMPLES}" \
            downsample_and_call "${tag}" "${target}" {} \
            < "${OUT_DIR}/metrics/samples_${tag}.txt"
        genotype_dataset "${tag}" "${target}"
    done

    # Every dataset has now been cut from the markdup BAMs.
    if ! ${KEEP_MARKDUP_BAM}; then
        log "  removing markdup BAMs (all datasets have consumed them)"
        rm -f "${OUT_DIR}"/markdup/*.markdup.bam "${OUT_DIR}"/markdup/*.markdup.bai
    fi
    log "PHASE 3 complete"
}

phase_downstream() {
    log "PHASE 4 — LD pruning, subsets and ROH inputs"
    for tag_target in "${DATASETS[@]}"; do
        local tag="${tag_target%%:*}"
        log "  dataset ${tag}"
        ld_prune_and_subset "${tag}"
        prepare_roh_input "${tag}"
    done
    log "PHASE 4 complete"
}

summary() {
    log "SUMMARY"
    for tag_target in "${DATASETS[@]}"; do
        # Assigned separately: within a single `local`, all words are expanded
        # before any assignment happens, so dd could not reference tag here.
        local tag target dd
        tag="${tag_target%%:*}"; target="${tag_target##*:}"; dd="${OUT_DIR}/${tag}"
        info "  ── ${tag} (${target}x, $(wc -l < "${OUT_DIR}/metrics/samples_${tag}.txt") samples) ──"
        for f in "${dd}/filtered/cohort.snps.final.vcf.gz" \
                 "${dd}/downstream/ld_pruned.snps.vcf.gz" \
                 "${dd}/downstream/ld_pruned_50K.snps.vcf.gz" \
                 "${dd}/downstream/ld_pruned_5K.snps.vcf.gz" \
                 "${dd}/roh/for_roh_dense_maf005.vcf.gz" \
                 "${dd}/roh/for_roh_dense_maf010.vcf.gz"; do
            if [[ -f "${f}" ]]; then
                info "    $(printf '%-46s %10s SNPs' "$(basename "${f}")" "$(bcftools index -n "${f}" 2>/dev/null || echo '?')")"
            else
                warn "    $(basename "${f}") MISSING"
            fi
        done
    done
    info "  disk $(disk_pct)% used, $(df -Ph "${WORK_DIR}" | awk 'NR==2{print $4}') free"
    log "Pipeline finished"
}

main() {
    preflight
    case "${PHASE}" in
        all)        phase_align; phase_assign; phase_call; phase_downstream; summary ;;
        align)      phase_align ;;
        assign)     phase_assign ;;
        call)       phase_assign; phase_call ;;
        downstream) phase_downstream; summary ;;
        *)          error "unknown phase '${PHASE}' (use: all|align|assign|call|downstream)" ;;
    esac
}

main
