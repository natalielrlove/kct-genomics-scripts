#!/bin/bash
# ============================================================
# 04b_mapq_filter.sh — Filter low-MAPQ reads from sorted BAMs
# ============================================================
# Filters out reads with mapping quality (MAPQ) < 30 from the
# sorted BAMs produced by 04_bwamem.sh, before read group
# assignment and duplicate marking in script 05.
#
# Applying this filter at the read level (before HaplotypeCaller)
# is more aggressive than post-hoc INFO/MQ filtering: reads are
# removed entirely so the variant caller never sees them.
# This is the approach recommended by the Kreiner Lab (ECEV35300
# Lab 2) for non-model organisms.
#
# MAPQ = 30 means 99.9% confidence the read is correctly placed.
# Reads below this threshold are likely multi-mappers from
# repetitive or poorly-assembled regions.
#
# Input:  sorted BAMs from 04_bwamem.sh
#           ~/kct_genomics/output_files/bwa_mem2_hap1/bam/*.hap1.sorted.bam
# Output: MAPQ-filtered BAMs
#           ~/kct_genomics/output_files/bwa_mem2_hap1/bam_mq30/*.hap1.sorted.mq30.bam
#
# These filtered BAMs replace the sorted BAMs as input to
# script 05 (AddOrReplaceReadGroups + MarkDuplicates).
#
# Parallelism:
#   - 10 samples processed simultaneously
#   - Each job uses 4 threads for samtools view + index
#   - Total CPU: ~40 cores
#
# Resume behavior:
#   - If a filtered BAM already exists for a sample, it is
#     skipped. Delete bam_mq30/ to force a full re-run.
#
# Requirements: conda environment "mapping" must be active
#   (contains samtools)
#   source /home/nlove/miniconda3/etc/profile.d/conda.sh
#   conda activate mapping
# ============================================================

set -euo pipefail

# --- Directory paths -----------------------------------------
IN="/home/nlove/kct_genomics/output_files/bwa_mem2_hap1/bam"
OUT="/home/nlove/kct_genomics/output_files/bwa_mem2_hap1/bam_mq30"

# --- Settings ------------------------------------------------
threads=4
parallel_jobs=10
mapq_threshold=30

# --- Create output directory ---------------------------------
mkdir -p "$OUT"

# --- Define per-sample function ------------------------------
process_sample() {
  local bam="$1"
  local out="$2"
  local threads="$3"
  local mapq="$4"

  local sample
  sample=$(basename "$bam" .hap1.sorted.bam)

  # Skip Undetermined
  if [[ "$sample" == Undetermined* ]]; then
    echo "Skipping $sample — Undetermined reads"
    return 0
  fi

  # Resume check
  if [[ -f "$out/${sample}.hap1.sorted.mq30.bam" ]]; then
    echo "Skipping $sample — filtered BAM already exists"
    return 0
  fi

  echo "Filtering: $sample"

  # Filter reads with MAPQ < threshold, output compressed BAM
  samtools view \
    -q "$mapq" \
    -b \
    -@ "$threads" \
    -o "$out/${sample}.hap1.sorted.mq30.bam" \
    "$bam"

  # Index the filtered BAM
  samtools index -@ "$threads" "$out/${sample}.hap1.sorted.mq30.bam"

  echo "Done: $sample"
}

export -f process_sample

# --- Run in parallel -----------------------------------------
echo "Starting MAPQ >= ${mapq_threshold} filtering for $(ls "$IN"/*.hap1.sorted.bam | wc -l) samples..."
echo "Running $parallel_jobs samples in parallel, $threads threads each."

parallel -j "$parallel_jobs" \
  process_sample {} "$OUT" "$threads" "$mapq_threshold" \
  ::: "$IN"/*.hap1.sorted.bam

echo ""
echo "All samples complete. Filtered BAMs are in: $OUT"
echo "Update script 05 input path to: $OUT"
