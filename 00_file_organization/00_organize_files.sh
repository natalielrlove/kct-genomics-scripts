#!/usr/bin/env bash
# =============================================================================
# 00_organize_files.sh
# Reorganize ~/kct_genomics/ to match the target architecture:
#
#   ~/kct_genomics/
#   ├── input_files/
#   │   ├── Gymnocladus_dioicus_male_assembly/
#   │   ├── NovaSeq/01.RawData/all_fastq/
#   │   └── reference/
#   ├── output_files/
#   └── scripts/
#
# Run from the server after pulling this repository.
# Review the DRY RUN output carefully before setting DRY_RUN=false.
# =============================================================================

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────
BASE="${HOME}/kct_genomics"

# Set to true to preview actions without making any changes.
# Set to false to execute.
DRY_RUN=false

# ── Helpers ───────────────────────────────────────────────────────────────────
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
drylog() {
    if [[ "${DRY_RUN}" == true ]]; then
        echo "[DRY RUN] $*"
    fi
}

run() {
    if [[ "${DRY_RUN}" == true ]]; then
        echo "[DRY RUN] would run: $*"
    else
        log "Running: $*"
        eval "$@"
    fi
}

# ── Pre-flight checks ─────────────────────────────────────────────────────────
log "Starting file reorganization for ${BASE}"
[[ "${DRY_RUN}" == true ]] && log "*** DRY RUN MODE — no files will be moved or deleted ***"

if [[ ! -d "${BASE}" ]]; then
    echo "ERROR: Base directory ${BASE} not found. Exiting."
    exit 1
fi

# ── Step 1: Create required directories ───────────────────────────────────────
log "Step 1: Creating directory structure"
run "mkdir -p '${BASE}/input_files/Gymnocladus_dioicus_male_assembly'"
run "mkdir -p '${BASE}/input_files/NovaSeq/01.RawData/all_fastq'"
run "mkdir -p '${BASE}/input_files/reference'"
run "mkdir -p '${BASE}/output_files'"
run "mkdir -p '${BASE}/scripts'"

# ── Step 2: Move Gymnocladus assembly into input_files/ ───────────────────────
log "Step 2: Moving Gymnocladus_dioicus_male_assembly into input_files/"

SRC_ASSEMBLY="${BASE}/Gymnocladus_dioicus_male_assembly"
DST_ASSEMBLY="${BASE}/input_files/Gymnocladus_dioicus_male_assembly"

if [[ -d "${SRC_ASSEMBLY}" ]]; then
    # Move contents (preserving subdirectories like plastid_assembly/ and mito_assembly/)
    run "rsync -av --progress '${SRC_ASSEMBLY}/' '${DST_ASSEMBLY}/'"
    if [[ "${DRY_RUN}" == false ]]; then
        log "Verifying move completed before removing source..."
        # Only remove source if destination has the same file count
        SRC_COUNT=$(find "${SRC_ASSEMBLY}" -type f | wc -l)
        DST_COUNT=$(find "${DST_ASSEMBLY}" -type f | wc -l)
        if [[ "${SRC_COUNT}" -eq "${DST_COUNT}" ]]; then
            log "File counts match (${SRC_COUNT}). Removing source directory."
            rm -rf "${SRC_ASSEMBLY}"
        else
            log "WARNING: Source has ${SRC_COUNT} files, destination has ${DST_COUNT}. Source NOT removed — check manually."
        fi
    fi
else
    log "WARNING: ${SRC_ASSEMBLY} not found — already moved, or path is wrong."
fi

# ── Step 3: Complete the NovaSeq fastq copy into input_files/ ─────────────────
log "Step 3: Syncing all .fq.gz files into input_files/NovaSeq/01.RawData/all_fastq/"

SRC_FASTQ="${BASE}/NovaSeq/01.RawData/all_fastq"
DST_FASTQ="${BASE}/input_files/NovaSeq/01.RawData/all_fastq"

if [[ -d "${SRC_FASTQ}" ]]; then
    SRC_FQ_COUNT=$(find "${SRC_FASTQ}" -name "*.fq.gz" | wc -l)
    DST_FQ_COUNT=$(find "${DST_FASTQ}" -name "*.fq.gz" 2>/dev/null | wc -l)
    log "  Source: ${SRC_FQ_COUNT} .fq.gz files"
    log "  Destination currently: ${DST_FQ_COUNT} .fq.gz files"

    # rsync --ignore-existing skips files already present; checksum validates integrity
    run "rsync -av --progress --checksum '${SRC_FASTQ}/' '${DST_FASTQ}/'"
else
    log "WARNING: ${SRC_FASTQ} not found. Skipping fastq sync."
fi

# ── Step 4: Remove duplicate top-level reference/ directory ───────────────────
log "Step 4: Checking for duplicate top-level reference/ directory"

TOP_REF="${BASE}/reference"
INPUT_REF="${BASE}/input_files/reference"

if [[ -d "${TOP_REF}" ]] && [[ -d "${INPUT_REF}" ]]; then
    TOP_COUNT=$(find "${TOP_REF}" -type f | wc -l)
    INPUT_COUNT=$(find "${INPUT_REF}" -type f | wc -l)
    log "  Top-level reference/ has ${TOP_COUNT} files"
    log "  input_files/reference/ has ${INPUT_COUNT} files"

    if [[ "${TOP_COUNT}" -eq "${INPUT_COUNT}" ]]; then
        log "  Counts match — safe to remove duplicate top-level reference/"
        run "rm -rf '${TOP_REF}'"
    else
        log "  WARNING: Counts differ. Top-level reference/ NOT removed — check manually."
    fi
elif [[ -d "${TOP_REF}" ]] && [[ ! -d "${INPUT_REF}" ]]; then
    log "  input_files/reference/ missing — moving top-level reference/ into input_files/"
    run "mv '${TOP_REF}' '${INPUT_REF}'"
else
    log "  No top-level reference/ found — nothing to do."
fi

# ── Step 5: Verify final structure ────────────────────────────────────────────
log "Step 5: Final directory structure"
if [[ "${DRY_RUN}" == false ]]; then
    find "${BASE}" -maxdepth 3 -not -path "*/\.*" | sort
else
    drylog "Would display final directory structure here."
fi

# ── Step 6: Post-copy cleanup (manual confirmation recommended) ───────────────
log ""
log "============================================================"
log "IMPORTANT: After verifying the final structure above,"
log "manually remove the original top-level NovaSeq/ directory:"
log ""
log "  rm -rf ${BASE}/NovaSeq/"
log ""
log "Do NOT automate this step — verify file counts first:"
log "  find ${BASE}/input_files/NovaSeq -name '*.fq.gz' | wc -l"
log "  find ${BASE}/NovaSeq -name '*.fq.gz' | wc -l"
log "============================================================"

log "Done."
