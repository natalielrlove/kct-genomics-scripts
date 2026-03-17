#!/bin/bash
# ============================================================
# 05_gatk_markdups.sh — Add read groups and mark PCR duplicates
# ============================================================
# This script prepares BAM files for variant calling. It runs
# two GATK tools on each sample, in order:
#
#   STEP 1 — AddOrReplaceReadGroups
#   --------------------------------
#   BAM files coming out of bwa-mem2 (script 04) contain the
#   aligned reads but no metadata about where or when they were
#   sequenced. GATK requires this metadata — called "read groups"
#   — to be embedded in the BAM before it will run HaplotypeCaller.
#
#   Read group fields we assign (parsed from the filename):
#     RGID (ID)      : unique identifier for this read group
#                      → set to sample_flowcell_lane (e.g. Gymno_2013_110_22V7CCLT4_L2)
#     RGSM (sample)  : biological sample name (e.g. Gymno_2013_110)
#     RGLB (library) : library prep barcode (e.g. CKDL250009714-1A)
#     RGPL (platform): always ILLUMINA for our data
#     RGPU (unit)    : flowcell + lane (e.g. 22V7CCLT4_L2)
#
#   STEP 2 — MarkDuplicates
#   -----------------------
#   During library preparation, PCR amplification can create
#   many identical copies of the same original DNA molecule.
#   If all these copies are counted as independent evidence
#   for a variant, variant calls become biased (inflated
#   confidence in common alleles, missed rare variants).
#
#   MarkDuplicates identifies reads that appear to be PCR
#   copies (same start position, same sequence, same strand)
#   and flags them in the BAM with a duplicate bit. They are
#   NOT deleted — just marked so that GATK can ignore them
#   during variant calling.
#
#   A metrics file is also written showing the duplication
#   rate for each sample. High duplication (>30-40%) may
#   indicate a library quality issue worth noting.
#
# Output per sample:
#   {sample}.rg.bam      : BAM with read groups added (intermediate)
#   {sample}.md.bam      : BAM with duplicates marked → input for script 06
#   {sample}.md.bai      : index for the duplicate-marked BAM
#   {sample}.metrics.txt : duplication statistics (check these!)
#   {sample}.rg.log      : GATK log for AddOrReplaceReadGroups
#   {sample}.markdup.log : GATK log for MarkDuplicates
#
# Parallelism:
#   - 10 samples processed simultaneously (parallel_jobs=10)
#   - MarkDuplicates is mostly single-threaded but memory-heavy
#   - Each job uses up to 8 GB RAM → 10 jobs = ~80 GB RAM total
#   - Total CPU: ~10 cores (one per job, GATK does not multithread here)
#
# Resume behavior:
#   - If a .md.bam already exists for a sample, that sample is
#     skipped automatically. Safe to re-run after interruption.
#
# Requirements: conda environment "gatk" must be active
#   (contains gatk4 and samtools)
# ============================================================

# Exit immediately if any command fails, if an unset variable
# is used, or if any command in a pipe fails
set -euo pipefail

# --- Activate conda environment ------------------------------
# The "gatk" environment contains GATK4 and samtools
source /home/nlove/miniconda3/etc/profile.d/conda.sh
conda activate gatk

# --- Directory paths -----------------------------------------
# Input: sorted BAM files produced by 04_bwamem.sh
IN="/home/nlove/kct_genomics/output_files/bwa_mem2_hap1/bam"

# Output: all GATK preprocessing results go here
OUT="/home/nlove/kct_genomics/output_files/gatk_preproc"

# --- Parallelism settings ------------------------------------
parallel_jobs=10  # Number of samples to process at the same time.
                  # MarkDuplicates is memory-intensive (~8 GB/job).
                  # 10 jobs × 8 GB = ~80 GB RAM. Increase if you have
                  # headroom, decrease if jobs crash with memory errors.

# --- Create output directories if they don't already exist ---
# rg/      = BAMs after read group assignment (intermediate step)
#            These can be deleted once md/ BAMs are confirmed good.
# md/      = BAMs after duplicate marking → used by HaplotypeCaller
# metrics/ = Per-sample duplication statistics (% duplicates, etc.)
# logs/    = GATK log output for each sample (useful for debugging)
mkdir -p "$OUT"/{rg,md,metrics,logs}

# --- Define the function that processes one sample -----------
# GNU parallel will call this function once per BAM file.
# It must be exported so parallel can access it in subprocesses.
process_sample() {

  # Receive the BAM filepath as the first argument from parallel
  local bam="$1"
  local out="$2"

  # Strip the path and ".hap1.sorted.bam" suffix to get the base name
  # e.g.: /path/to/Gymno_2013_110_CKDL250009714-1A_22V7CCLT4_L2.hap1.sorted.bam
  #     →                Gymno_2013_110_CKDL250009714-1A_22V7CCLT4_L2
  local base
  base=$(basename "$bam" .hap1.sorted.bam)

  # --- Undetermined skip check --------------------------------
  # Undetermined BAMs come from reads that couldn't be assigned
  # to any sample during demultiplexing. Skip them.
  if [[ "$base" == Undetermined* ]]; then
    echo "Skipping $base — Undetermined reads, not a real sample"
    return 0
  fi

  # --- Resume/skip check -------------------------------------
  # If a duplicate-marked BAM already exists, skip this sample.
  # This lets you safely re-run after an interruption without
  # re-processing samples that already finished.
  if [[ -f "$out/md/${base}.md.bam" ]]; then
    echo "Skipping $base — md.bam already exists"
    return 0
  fi

  # --- Parse read group fields from the filename --------------
  # Filenames follow the pattern:
  #   Gymno_2013_110_CKDL250009714-1A_22V7CCLT4_L2.hap1.sorted.bam
  #   |---- SM ----|  |---- LB ----|  |--- PU ---|
  #
  # SM = sample name = first 3 underscore-delimited fields
  #   cut -d_ : split on underscore
  #   -f1-3   : take fields 1 through 3
  #   result  : "Gymno_2013_110"
  local sm
  sm=$(echo "$base" | cut -d_ -f1-3)

  # LB = library prep barcode = 4th underscore-delimited field
  #   result: "CKDL250009714-1A"
  local lb
  lb=$(echo "$base" | cut -d_ -f4)

  # PU = platform unit = last two underscore-delimited fields
  #   awk -F_ splits on underscore; NF is total number of fields
  #   $(NF-1) = second-to-last field (flowcell ID)
  #   $NF     = last field (lane)
  #   result  : "22V7CCLT4_L2"
  local pu
  pu=$(echo "$base" | awk -F_ '{print $(NF-1)"_"$NF}')

  # RGID = read group ID = sample name + platform unit
  #   Must be unique across all samples in a multi-sample project
  #   result: "Gymno_2013_110_22V7CCLT4_L2"
  local rgid="${sm}_${pu}"

  echo "Processing: $base"
  echo "  SM=$sm  LB=$lb  PU=$pu  RGID=$rgid"

  # --- Step 1: Add Read Groups --------------------------------
  # Embeds SM, LB, PL, PU, RGID as tags in every read in the BAM.
  # GATK HaplotypeCaller will refuse to run without these tags.
  #
  # --java-options "-Xmx8g" : tell Java it can use up to 8 GB RAM
  # -I : input BAM (sorted, no read groups)
  # -O : output BAM (same reads, with RG tags embedded)
  # -RGID / -RGLB / -RGPL / -RGPU / -RGSM : the tag values to add
  # 2> : redirect GATK's progress messages to a log file
  gatk --java-options "-Xmx8g" AddOrReplaceReadGroups \
    -I "$bam" \
    -O "$out/rg/${base}.rg.bam" \
    -RGID "$rgid" \
    -RGLB "$lb" \
    -RGPL ILLUMINA \
    -RGPU "$pu" \
    -RGSM "$sm" \
    2> "$out/logs/${base}.rg.log"

  # --- Step 2: Mark Duplicates --------------------------------
  # Scans the read-group-tagged BAM and identifies reads that
  # appear to be PCR duplicates of the same original molecule
  # (same chromosome, same start position, same orientation).
  # Duplicates are flagged with a bit in the BAM flags — NOT
  # deleted. GATK will skip flagged duplicates automatically.
  #
  # --java-options "-Xmx8g" : up to 8 GB RAM for Java
  # -I : input BAM (with read groups from step 1)
  # -O : output BAM (duplicates flagged)
  # -M : metrics file — reports % duplication rate per sample.
  #      Check these after the run: >30-40% duplication may
  #      indicate a library prep issue worth noting in your methods.
  # --CREATE_INDEX true : automatically creates a .bai index file
  #                       alongside the output BAM (required by GATK)
  # 2> : redirect GATK's progress messages to a log file
  gatk --java-options "-Xmx8g" MarkDuplicates \
    -I "$out/rg/${base}.rg.bam" \
    -O "$out/md/${base}.md.bam" \
    -M "$out/metrics/${base}.metrics.txt" \
    --CREATE_INDEX true \
    2> "$out/logs/${base}.markdup.log"

  echo "Done: $base"
}

# Export the function so GNU parallel can call it in subprocesses
export -f process_sample

# --- Run in parallel -----------------------------------------
# GNU parallel distributes the BAM list across parallel_jobs slots.
# {} is replaced with each BAM filepath in turn.
echo "Starting read group assignment and duplicate marking..."
echo "Samples found: $(ls "$IN"/*.hap1.sorted.bam | wc -l)"
echo "Running $parallel_jobs samples in parallel."

parallel -j "$parallel_jobs" \
  process_sample {} "$OUT" \
  ::: "$IN"/*.hap1.sorted.bam

echo ""
echo "All samples complete."
echo "Duplicate-marked BAMs are in:    $OUT/md/"
echo "Duplication statistics are in:   $OUT/metrics/"
echo "Check metrics files for duplication rates before proceeding to script 06."
