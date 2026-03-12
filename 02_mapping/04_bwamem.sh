#!/bin/bash
# ============================================================
# 04_bwamem.sh — Align trimmed reads to reference genome
# ============================================================
# This script maps cleaned, trimmed reads (from fastp) to the
# Gymnocladus dioicus reference genome (haplotype 1) using
# bwa-mem2, a fast short-read aligner.
#
# For each sample it produces:
#   - A sorted BAM file: the aligned reads, sorted by position
#   - A BAM index (.bai): required to quickly access the BAM
#   - A flagstat file: alignment summary statistics per sample
#
# Parallelism:
#   - 2 samples are processed simultaneously (parallel_jobs=2)
#   - Each job uses 15 CPU threads for alignment (threads_bwa=15)
#   - Each job uses 2 CPU threads for sorting (threads_sort=2)
#   - Total CPU usage: 2 × (15 + 2) = 34 cores
#
# Resume behavior:
#   - If a sorted BAM already exists for a sample, that sample
#     is skipped automatically. This means if the job is
#     interrupted, you can re-run the script and it will pick
#     up where it left off without redoing finished samples.
#
# Requirements: conda environment "mapping" must be active
#   (contains bwa-mem2, samtools, and GNU parallel)
# ============================================================

# Exit immediately if any command fails, if an unset variable
# is used, or if any command in a pipe fails
set -euo pipefail

# --- Activate conda environment ------------------------------
# The "mapping" environment contains bwa-mem2, samtools, and parallel
source /home/nlove/miniconda3/etc/profile.d/conda.sh
conda activate mapping

# --- Directory paths -----------------------------------------
# Reference genome (haplotype 1) — already indexed, no need to re-index
REF="/home/nlove/kct_genomics/input_files/reference/Gymnocladus_dioicus_M_hap1.fa"

# Input: trimmed FASTQ files produced by 03_fastp.sh
IN="/home/nlove/kct_genomics/output_files/fastp/trimmed"

# Output: all alignment results go here
OUT="/home/nlove/kct_genomics/output_files/bwa_mem2_hap1"

# --- Thread and parallelism settings -------------------------
threads_bwa=15    # CPU threads used by bwa-mem2 per sample
threads_sort=2    # CPU threads used by samtools sort per sample
parallel_jobs=2   # Number of samples to process at the same time

# --- Create output directories if they don't already exist ---
# bam/   = sorted BAM alignment files and their indexes
# stats/ = per-sample flagstat alignment summary files
# logs/  = per-sample bwa-mem2 log files (useful for debugging)
mkdir -p "$OUT"/{bam,stats,logs}

# --- Load reference genome into shared memory ----------------
# bwa-mem2 normally loads the reference index from disk for every
# sample. With "shm", it loads it into RAM once and all parallel
# jobs share the same copy — this saves time and memory.
echo "Loading reference genome into shared memory..."
bwa-mem2 shm "$REF"

# --- Clean up shared memory when the script finishes ---------
# The "trap" command ensures shared memory is released when the
# script exits (whether it finishes normally or crashes).
# Without this, the reference would stay loaded in RAM indefinitely.
trap 'echo "Releasing shared memory..."; bwa-mem2 shm -d "$REF"' EXIT

# --- Define the function that processes one sample -----------
# GNU parallel will call this function once per sample.
# It must be defined as a function and exported so parallel
# can access it in its subprocesses.
process_sample() {

  # Receive the R1 filepath as the first argument from parallel
  local r1="$1"

  # Construct the R2 filename by replacing _1.trim.fq.gz with _2.trim.fq.gz
  # e.g. Gymno_2013_109_..._1.trim.fq.gz → Gymno_2013_109_..._2.trim.fq.gz
  local r2="${r1%_1.trim.fq.gz}_2.trim.fq.gz"

  # Extract the sample name (no path, no suffix)
  # e.g. /path/to/Gymno_2013_109_CKDL..._1.trim.fq.gz → Gymno_2013_109_CKDL...
  local sample
  sample=$(basename "$r1" _1.trim.fq.gz)

  # Receive remaining arguments passed in by parallel
  local out="$2"
  local ref="$3"
  local threads_bwa="$4"
  local threads_sort="$5"

  # --- Undetermined skip check --------------------------------
  # Undetermined files contain reads that couldn't be assigned to
  # any sample during demultiplexing — they are not real samples
  # and should not be aligned to the reference genome.
  if [[ "$sample" == Undetermined* ]]; then
    echo "Skipping $sample — Undetermined reads, not a real sample"
    return 0
  fi

  # --- Resume/skip check -------------------------------------
  # If a sorted BAM already exists for this sample, skip it.
  # This allows the script to be safely re-run after interruption
  # without re-aligning samples that already completed.
  if [[ -f "$out/bam/${sample}.hap1.sorted.bam" ]]; then
    echo "Skipping $sample — BAM already exists"
    return 0
  fi

  echo "Aligning: $sample"

  # --- Align reads and sort ----------------------------------
  # This is a pipeline (|) — the output of bwa-mem2 is piped
  # directly into samtools sort without writing a large
  # intermediate file to disk.
  #
  # bwa-mem2 mem: aligns paired-end reads to the reference
  #   -t $threads_bwa : use this many CPU threads
  #   "$ref"          : path to reference genome (loaded in shared memory)
  #   "$r1" "$r2"     : forward and reverse read files
  #   2> logs/...     : save bwa-mem2 progress messages to a log file
  #
  # samtools sort: sorts aligned reads by genome position
  #   -@ $threads_sort          : use this many extra threads
  #   -T $out/bam/${sample}.tmp : temporary file prefix for sorting
  #   -o $out/bam/...sorted.bam : output sorted BAM file
  #   -                         : read input from the pipe (stdout of bwa-mem2)
  bwa-mem2 mem -t "$threads_bwa" \
    "$ref" "$r1" "$r2" \
    2> "$out/logs/${sample}.bwa.log" \
    | samtools sort \
        -@ "$threads_sort" \
        -T "$out/bam/${sample}.tmp" \
        -o "$out/bam/${sample}.hap1.sorted.bam" \
        -

  # --- Index the sorted BAM ----------------------------------
  # Creates a .bai index file alongside the BAM.
  # This is required by downstream tools (e.g. GATK) to quickly
  # jump to any position in the genome without reading the whole file.
  samtools index -@ "$threads_bwa" "$out/bam/${sample}.hap1.sorted.bam"

  # --- Generate alignment statistics -------------------------
  # flagstat counts how many reads mapped, how many are duplicates,
  # how many are properly paired, etc. One .txt file per sample.
  samtools flagstat -@ "$threads_bwa" "$out/bam/${sample}.hap1.sorted.bam" \
    > "$out/stats/${sample}.hap1.flagstat.txt"

  echo "Done: $sample"
}

# Export the function so GNU parallel can call it in subprocesses
export -f process_sample

# --- Run alignment in parallel -------------------------------
# GNU parallel distributes the sample list across parallel_jobs slots.
# {} is replaced with each R1 filepath in turn.
# The extra arguments (out, ref, threads) are passed to each call.
echo "Starting alignment of $(ls "$IN"/*_1.trim.fq.gz | wc -l) samples..."
echo "Running $parallel_jobs samples in parallel, $threads_bwa threads each."

parallel -j "$parallel_jobs" \
  process_sample {} "$OUT" "$REF" "$threads_bwa" "$threads_sort" \
  ::: "$IN"/*_1.trim.fq.gz

echo "All samples complete. BAM files are in: $OUT/bam/"
