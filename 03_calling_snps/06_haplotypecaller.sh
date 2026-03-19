#!/bin/bash
# ============================================================
# 06_haplotypecaller.sh — Per-sample variant calling (gVCF mode)
# ============================================================
# This script runs GATK HaplotypeCaller on each sample to
# identify locations in the genome where reads differ from
# the reference — these differences are potential SNPs or
# small insertions/deletions (indels).
#
# WHY gVCF MODE?
# --------------
# HaplotypeCaller is run here in "gVCF mode" (-ERC GVCF).
# Instead of producing a final list of variants, it produces
# a per-sample gVCF file that records:
#   - Sites where a variant WAS detected (with confidence scores)
#   - Sites where NO variant was detected (with evidence that
#     the reference allele was actually observed there)
#
# This "evidence at every site" format is required for the
# next step: joint genotyping across all samples together
# (GenomicsDBImport + GenotypeGVCFs, scripts 07-08).
# Joint genotyping is more accurate than calling each sample
# independently because it uses information from all samples
# to distinguish true rare variants from sequencing errors.
#
# Think of it like this:
#   Script 06 (this script) = "what did each individual see?"
#   Script 07-08            = "what does the whole population agree on?"
#
# Output per sample:
#   {sample}.g.vcf.gz     : the per-sample gVCF file (compressed)
#   {sample}.g.vcf.gz.tbi : index for the gVCF (required for script 07)
#   {sample}.hc.log       : HaplotypeCaller log (useful for debugging)
#
# Runtime note:
#   HaplotypeCaller is the slowest step in this pipeline.
#   Expect ~8 hours per sample on a single thread. With all samples
#   running in parallel at 3 threads each, total wall time is
#   reduced to ~8 hours for all samples at once (using 30 cores).
#
# Parallelism:
#   - 10 samples processed simultaneously (parallel_jobs=10)
#   - Each job uses 3 CPU threads (threads_hc=3)
#   - Total CPU: 10 × 3 = 30 cores
#   - Each job uses ~16 GB RAM → 10 jobs = ~160 GB RAM total
#     (adjust parallel_jobs down if the server runs out of memory)
#
# Resume behavior:
#   - If a .g.vcf.gz already exists for a sample, that sample
#     is skipped. Safe to re-run after interruption.
#
# Requirements: conda environment "gatk" must be active
#   (contains gatk4)
# ============================================================

# Exit immediately if any command fails, if an unset variable
# is used, or if any command in a pipe fails
set -euo pipefail

# --- Activate conda environment ------------------------------
# The "gatk" environment contains GATK4.
# Activate manually before running this script:
#   source /home/nlove/miniconda3/etc/profile.d/conda.sh
#   conda activate gatk
# source /home/nlove/miniconda3/etc/profile.d/conda.sh
# conda activate gatk

# --- Directory paths -----------------------------------------
# Reference genome (haplotype 1) — same reference used in script 04
REF="/home/nlove/kct_genomics/input_files/reference/Gymnocladus_dioicus_M_hap1.fa"

# Input: duplicate-marked BAMs produced by 05_gatk_markdups.sh
IN="/home/nlove/kct_genomics/output_files/gatk_preproc/md"

# Output: per-sample gVCF files go here
OUT="/home/nlove/kct_genomics/output_files/gatk_gvcf_hap1"

# --- Thread and parallelism settings -------------------------
threads_hc=3      # CPU threads used by HaplotypeCaller per sample
                  # (--native-pair-hmm-threads controls this)
parallel_jobs=10  # Number of samples to process at the same time
                  # 10 jobs × 3 threads = 30 cores total
                  # 10 jobs × ~16 GB RAM = ~160 GB RAM total
                  # Reduce parallel_jobs if you hit memory limits

# --- Create output directories if they don't already exist ---
# gvcf/ = per-sample gVCF files (input for joint genotyping)
# logs/ = per-sample HaplotypeCaller log files
mkdir -p "$OUT"/{gvcf,logs}

# --- Define the function that processes one sample -----------
# GNU parallel will call this function once per BAM file.
# It must be exported so parallel can access it in subprocesses.
process_sample() {

  # Receive arguments passed in by parallel
  local bam="$1"
  local out="$2"
  local ref="$3"
  local threads_hc="$4"

  # Strip the path and ".md.bam" suffix to get the base sample name
  # e.g.: /path/to/Gymno_2013_110_CKDL250009714-1A_22V7CCLT4_L2.md.bam
  #     →           Gymno_2013_110_CKDL250009714-1A_22V7CCLT4_L2
  local base
  base=$(basename "$bam" .md.bam)

  # --- Undetermined skip check --------------------------------
  if [[ "$base" == Undetermined* ]]; then
    echo "Skipping $base — Undetermined reads, not a real sample"
    return 0
  fi

  # --- Resume/skip check -------------------------------------
  # If a gVCF already exists for this sample, skip it.
  # Allows safe re-runs after interruption.
  if [[ -f "$out/gvcf/${base}.g.vcf.gz" ]]; then
    echo "Skipping $base — gVCF already exists"
    return 0
  fi

  echo "Calling variants: $base"

  # --- Run HaplotypeCaller ------------------------------------
  # HaplotypeCaller uses a local de-novo assembly approach:
  # for each region of the genome where reads suggest a variant,
  # it reassembles the reads into candidate haplotypes and
  # compares them against the reference to call variants.
  #
  # --java-options "-Xmx16g"
  #   Allow Java to use up to 16 GB RAM. HaplotypeCaller is
  #   memory-intensive, especially for high-coverage (30×) data.
  #
  # -R : reference genome FASTA (must be the same one used in bwa-mem2)
  #
  # -I : input BAM (duplicate-marked, with read groups)
  #
  # -O : output gVCF file (compressed with bgzip, .gz extension)
  #      The .tbi index is created automatically alongside it.
  #
  # -ERC GVCF
  #   "Emit Reference Confidence" in GVCF mode.
  #   This is what produces the per-site confidence records
  #   needed for joint genotyping (see header comment above).
  #   Without this flag, you'd get a standard VCF with only
  #   variant sites — which cannot be used for joint calling.
  #
  # --native-pair-hmm-threads
  #   Number of threads for the pair-HMM likelihood calculation
  #   (the most compute-intensive part of HaplotypeCaller).
  #   This is the correct way to multithread HaplotypeCaller —
  #   NOT -nct, which is deprecated in GATK4.
  #
  # 2> : redirect GATK's verbose progress output to a log file
  #      (it prints a LOT — you don't want this flooding the screen)
  gatk --java-options "-Xmx16g" HaplotypeCaller \
    -R "$ref" \
    -I "$bam" \
    -O "$out/gvcf/${base}.g.vcf.gz" \
    -ERC GVCF \
    --native-pair-hmm-threads "$threads_hc" \
    2> "$out/logs/${base}.hc.log"

  echo "Done: $base → $out/gvcf/${base}.g.vcf.gz"
}

# Export the function so GNU parallel can call it in subprocesses
export -f process_sample

# --- Run HaplotypeCaller in parallel -------------------------
# GNU parallel distributes the BAM list across parallel_jobs slots.
# Each slot runs process_sample with its BAM + shared arguments.
echo "Starting HaplotypeCaller variant calling..."
echo "Samples found: $(ls "$IN"/*.md.bam | wc -l)"
echo "Running $parallel_jobs samples in parallel, $threads_hc threads each."
echo "Expected runtime: ~8 hours (parallel across all samples)."
echo ""

parallel -j "$parallel_jobs" \
  process_sample {} "$OUT" "$REF" "$threads_hc" \
  ::: "$IN"/*.md.bam

echo ""
echo "All samples complete."
echo "Per-sample gVCF files are in: $OUT/gvcf/"
echo "Next step: run 07_genomicsdb_import.sh to combine gVCFs for joint genotyping."
