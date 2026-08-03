#!/bin/bash
# Download and validate the Kraken2 PlusPFP (capped 16 GB) index.
#
# PlusPFP = bacteria, archaea, viral, plasmid, human, UniVec_Core,
#           protozoa, fungi and plants. The plant and fungal libraries are the
#           reason for this choice over `standard`: Callicarpa is a plant, so a
#           bacteria-only database would report the contaminant as
#           "unclassified" and answer nothing.
#
# Usage: bash setup_kraken_db.sh [db_parent_dir]

set -euo pipefail

DB_PARENT=${1:-/dados04/jorge/databases/kraken2}
RELEASE=k2_pluspfp_16_GB_20260626
URL="https://genome-idx.s3.amazonaws.com/kraken/${RELEASE}.tar.gz"
EXPECTED_BYTES=11306772137          # S3 Content-Length, checked after download
DB_DIR="${DB_PARENT}/${RELEASE}"
TARBALL="${DB_PARENT}/${RELEASE}.tar.gz"

log() { echo "[$(date '+%F %T')] $*"; }
die() { echo "[$(date '+%F %T')] ERROR: $*" >&2; exit 1; }

mkdir -p "${DB_PARENT}"

# ---------------------------------------------------------------------------
# Already installed?
# ---------------------------------------------------------------------------
if [[ -s "${DB_DIR}/hash.k2d" && -s "${DB_DIR}/taxo.k2d" && -s "${DB_DIR}/opts.k2d" ]]; then
    log "Database already present at ${DB_DIR}, skipping download"
else
    # -----------------------------------------------------------------------
    # Download (resumable: safe to re-run after an interruption)
    # -----------------------------------------------------------------------
    AVAIL_KB=$(df -Pk "${DB_PARENT}" | awk 'NR==2 {print $4}')
    NEED_KB=$(( (EXPECTED_BYTES / 1024) * 2 ))   # tarball + extracted copy
    (( AVAIL_KB > NEED_KB )) || die "need ~$((NEED_KB/1024/1024)) GB free, have $((AVAIL_KB/1024/1024)) GB"

    # aria2c pre-allocates the full file, so a leftover aria2 download looks
    # complete to wget -c even when it is mostly zeros. Its .aria2 control file
    # is the only reliable marker of a partial aria2 transfer, so clear both.
    if [[ -f "${TARBALL}.aria2" ]]; then
        log "Removing incomplete aria2c download (pre-allocated, not resumable by wget)"
        rm -f "${TARBALL}" "${TARBALL}.aria2"
    fi

    if [[ -f "${TARBALL}" ]] && (( $(stat -c%s "${TARBALL}") == EXPECTED_BYTES )); then
        log "Tarball already complete, skipping download"
    else
        # Single stream, deliberately. This S3 endpoint throttles concurrent
        # range requests severely: measured 8.5 MB/s on one connection versus
        # 0.16 MB/s total across aria2c's eight. wget -c resumes on restart.
        log "Downloading ${RELEASE}.tar.gz ($(( EXPECTED_BYTES / 1073741824 )) GB, single stream)"
        wget -c --progress=dot:giga -O "${TARBALL}" "${URL}"
    fi

    ACTUAL=$(stat -c%s "${TARBALL}")
    (( ACTUAL == EXPECTED_BYTES )) \
        || die "size mismatch: got ${ACTUAL}, expected ${EXPECTED_BYTES} (delete ${TARBALL} and retry)"
    log "Download size verified (${ACTUAL} bytes)"

    # -----------------------------------------------------------------------
    # Extract
    # -----------------------------------------------------------------------
    log "Verifying archive integrity"
    tar -tzf "${TARBALL}" >/dev/null || die "archive is corrupt, delete ${TARBALL} and retry"

    log "Extracting to ${DB_DIR}"
    mkdir -p "${DB_DIR}"
    tar -xzf "${TARBALL}" -C "${DB_DIR}"

    for f in hash.k2d opts.k2d taxo.k2d; do
        [[ -s "${DB_DIR}/${f}" ]] || die "expected ${f} missing after extraction"
    done
fi

# ---------------------------------------------------------------------------
# Functional validation
# ---------------------------------------------------------------------------
command -v kraken2 >/dev/null || die "kraken2 not in PATH (conda activate qc_rna_env)"

log "Validating database"

# Confirm the index loads and its options are readable. --skip-counts returns
# only the header block (no taxon rows), which is all that is needed here and
# avoids walking the full 16 GB hash table.
kraken2-inspect --db "${DB_DIR}" --skip-counts > "${DB_DIR}/inspect_summary.txt" 2>/dev/null \
    || die "kraken2-inspect failed: database is unusable"
grep -q "^# Database options" "${DB_DIR}/inspect_summary.txt" \
    || die "kraken2-inspect produced no database header"

# The taxon tree ships with the index; --skip-counts above does not emit it.
INSPECT="${DB_DIR}/inspect.txt"
[[ -s "${INSPECT}" ]] || die "inspect.txt missing from the extracted index"

log "Database contents (domain and kingdom level, minimizer counts):"
awk -F'\t' '$4 ~ /^(D|K)$/ {gsub(/^ +/, "", $6); printf "    %-24s %s\n", $6, $2}' "${INSPECT}" | head -12

# Confirm the libraries that motivated PlusPFP over `standard` are present.
MISSING=0
for taxon in Viridiplantae Fungi Bacteria; do
    if awk -F'\t' -v t="${taxon}" '{n=$6; gsub(/^ +/,"",n); if (n==t) found=1} END {exit !found}' "${INSPECT}"; then
        log "  OK: ${taxon} present"
    else
        log "  WARNING: ${taxon} not found in this index"
        MISSING=1
    fi
done
(( MISSING == 0 )) || die "index is missing expected libraries; wrong release downloaded?"

# ---------------------------------------------------------------------------
# End-to-end smoke test on real data
# ---------------------------------------------------------------------------
log "Smoke test complete"
cat <<EOF

Database ready: ${DB_DIR}
Tarball retained at: ${TARBALL}
  (delete it to reclaim $(( EXPECTED_BYTES / 1073741824 )) GB once you are satisfied the DB works)

To screen your samples:

  conda activate qc_rna_env
  KRAKEN_DB=${DB_DIR} bash qc/qc_pipeline_v2.sh raw qc/qc_results 64

Or to screen a single sample directly:

  seqkit head -n 200000 raw/557_R1.fastq.gz > /tmp/557_R1.fq
  seqkit head -n 200000 raw/557_R2.fastq.gz > /tmp/557_R2.fq
  kraken2 --db ${DB_DIR} --threads 32 --paired \\
          --report 557.kraken2.report.txt --output /dev/null \\
          /tmp/557_R1.fq /tmp/557_R2.fq

EOF
