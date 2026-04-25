#!/bin/bash
# ============================================================
# 04_variantcall.sh — Align reads to supercontig reference and call variants
# ============================================================
# Per-sample pipeline: BWA align → AddOrReplaceReadGroups →
# MarkDuplicates → HaplotypeCaller (GVCF mode).
#
# This mirrors the genome pipeline (scripts 04 + 05 + 06) but
# maps reads to the clean supercontig reference (~1.9 Mb across
# 928 loci) rather than the full Gymnocladus genome (710 Mb).
# Combining the three steps into one script follows the Slimp
# et al. (2021) variantcall.sh structure.
#
# KEY DIFFERENCE FROM THE GENOME PIPELINE
# ----------------------------------------
# Only ~0.6–0.9% of WGS reads map to the supercontig reference
# (the same on-target fraction observed in HybPiper assembly).
# BWA still reads all raw FASTQ data to find these reads, so
# alignment runtime is similar to the genome pipeline. However:
#   - Sorted BAMs are tiny (~100–200 MB vs ~15 GB for the genome)
#   - MarkDuplicates runs in seconds (small BAM)
#   - HaplotypeCaller is dramatically faster (~minutes vs ~17h)
#     because the reference is 375× smaller
#
# READ GROUPS (matching genome pipeline script 05)
# -------------------------------------------------
# Read groups are added with GATK AddOrReplaceReadGroups after
# alignment, parsing SM/LB/PU from the filename — the same
# approach as script 05_gatk_markdups.sh. We do NOT use
# MergeBamAlignment (as in Slimp et al.) because our raw input
# is already full-genome WGS data; creating a uBAM for all
# reads would generate ~15–20 GB intermediate files per sample
# (~1 TB total) for only ~0.6–0.9% on-target reads.
#
# OUTPUT PER SAMPLE
# -----------------
#   bam/{sample}.sc.sorted.bam         : aligned, sorted BAM
#   bam/{sample}.sc.sorted.bam.bai     : BAM index
#   stats/{sample}.sc.flagstat.txt     : alignment summary
#   rg/{sample}.sc.rg.bam              : BAM with read groups (intermediate)
#   md/{sample}.sc.md.bam              : duplicate-marked BAM
#   md/{sample}.sc.md.bai              : index for duplicate-marked BAM
#   metrics/{sample}.sc.metrics.txt    : duplication rate statistics
#   gvcf/{sample}.sc.g.vcf.gz          : per-sample gVCF
#   gvcf/{sample}.sc.g.vcf.gz.tbi      : gVCF index
#   logs/                              : per-sample log files
#
# PARALLELISM
# -----------
# GNU parallel runs N_JOBS samples simultaneously.
# Each job uses THREADS_BWA threads for alignment.
# HaplotypeCaller uses THREADS_HC threads per sample.
# Total CPU = N_JOBS × THREADS_BWA (BWA step, rate-limiting).
#
# RUNTIME ESTIMATE
# ----------------
# BWA: ~30–60 min per sample (dominated by reading raw FASTQ)
# MarkDuplicates: ~1–2 min per sample (small BAM)
# HaplotypeCaller: ~5–15 min per sample (small reference)
# With N_JOBS=4 → ~2 hours total wall time
#
# Requirements: conda environment "gatk" must be active
#   (must contain bwa, samtools, gatk4, and GNU parallel)
#   source /home/nlove/miniconda3/etc/profile.d/conda.sh
#   conda activate gatk
# ============================================================

set -euo pipefail

# ============================================================
# --- Configuration -------------------------------------------
# ============================================================

# Clean supercontig reference built by 03_build_supercontig_ref.sh
REF="/home/nlove/kct_genomics/output_files/hybpiper/gymno_supercontigs_clean.fasta"

# Input: trimmed FASTQ files (same as genome pipeline)
IN="/home/nlove/kct_genomics/output_files/fastp/trimmed"

# Output directory for all variant calling results
OUT="/home/nlove/kct_genomics/output_files/hybpiper_variantcall"

# --- Thread and parallelism settings -------------------------
N_JOBS=4          # Samples to process simultaneously.
                  # BWA reads all raw FASTQ per sample — each job
                  # is I/O-heavy. 4 jobs × 8 threads = 32 cores.
THREADS_BWA=8     # BWA threads per sample
THREADS_HC=3      # HaplotypeCaller threads per sample (pair-HMM)
                  # Same as genome pipeline script 06.

# ============================================================
# --- Setup ---------------------------------------------------
# ============================================================

mkdir -p "$OUT"/{bam,rg,md,gvcf,metrics,stats,logs}

# ============================================================
# --- Per-sample function -------------------------------------
# ============================================================

process_sample() {

    local r1="$1"
    local in="$2"
    local out="$3"
    local ref="$4"
    local threads_bwa="$5"
    local threads_hc="$6"

    local r2="${r1/_1.trim.fq.gz/_2.trim.fq.gz}"

    # Extract full base name, then clean sample name (Gymno_YYYY_NNN)
    local base
    base=$(basename "$r1" _1.trim.fq.gz)
    local sample
    sample=$(echo "$base" | cut -d'_' -f1-3)

    # --- Skip Undetermined reads --------------------------------
    if [[ "$base" == Undetermined* ]]; then
        echo "Skipping $base — Undetermined reads"
        return 0
    fi

    # --- Resume/skip check -------------------------------------
    # Skip if the final gVCF already exists for this sample
    if [[ -f "$out/gvcf/${sample}.sc.g.vcf.gz" ]]; then
        echo "Skipping $sample — gVCF already exists"
        return 0
    fi

    # --- Parse read group fields from filename -----------------
    # Follows the same logic as genome pipeline script 05.
    # Filename: Gymno_2013_110_CKDL250009714-1A_22V7CCLT4_L2_1.trim.fq.gz
    #   SM  = Gymno_2013_110
    #   LB  = CKDL250009714-1A
    #   PU  = 22V7CCLT4_L2
    #   RGID = Gymno_2013_110_22V7CCLT4_L2
    local sm="$sample"
    local lb
    lb=$(echo "$base" | cut -d'_' -f4)
    local pu
    pu=$(echo "$base" | awk -F_ '{print $(NF-1)"_"$NF}')
    local rgid="${sm}_${pu}"

    echo "========================================"
    echo "Processing: $sample"
    echo "  SM=$sm  LB=$lb  PU=$pu  RGID=$rgid"
    echo "========================================"

    # ===========================================================
    # STEP A: Align with BWA
    # ===========================================================
    # Maps trimmed reads to the clean supercontig reference.
    # Pipes directly to samtools sort to avoid a large temp file.
    # Uses bwa mem (not bwa-mem2) because the supercontig
    # reference was indexed with bwa in 03_build_supercontig_ref.sh.

    echo "[$sample] Aligning with BWA..."

    bwa mem -t "$threads_bwa" \
        "$ref" "$r1" "$r2" \
        2> "$out/logs/${sample}.bwa.log" \
    | samtools sort \
        -@ 2 \
        -T "$out/bam/${sample}.tmp" \
        -o "$out/bam/${sample}.sc.sorted.bam" \
        -

    samtools index "$out/bam/${sample}.sc.sorted.bam"

    samtools flagstat "$out/bam/${sample}.sc.sorted.bam" \
        > "$out/stats/${sample}.sc.flagstat.txt"

    echo "[$sample] BWA done. Flagstat written."

    # ===========================================================
    # STEP B: Add Read Groups
    # ===========================================================
    # GATK HaplotypeCaller requires read group tags (SM, LB, PL,
    # PU, RGID) embedded in the BAM. Parsed from the filename,
    # matching genome pipeline script 05_gatk_markdups.sh.

    echo "[$sample] Adding read groups..."

    gatk --java-options "-Xmx8g" AddOrReplaceReadGroups \
        -I "$out/bam/${sample}.sc.sorted.bam" \
        -O "$out/rg/${sample}.sc.rg.bam" \
        -RGID "$rgid" \
        -RGLB "$lb" \
        -RGPL ILLUMINA \
        -RGPU "$pu" \
        -RGSM "$sm" \
        2> "$out/logs/${sample}.rg.log"

    # ===========================================================
    # STEP C: Mark Duplicates
    # ===========================================================
    # Flags PCR duplicate reads by position. Duplicates are marked
    # not deleted — GATK ignores them automatically.
    # Check metrics files: duplication rates for WGS-based virtual
    # capture may differ from true capture data.

    echo "[$sample] Marking duplicates..."

    gatk --java-options "-Xmx8g" MarkDuplicates \
        -I "$out/rg/${sample}.sc.rg.bam" \
        -O "$out/md/${sample}.sc.md.bam" \
        -M "$out/metrics/${sample}.sc.metrics.txt" \
        --CREATE_INDEX true \
        2> "$out/logs/${sample}.markdup.log"

    # ===========================================================
    # STEP D: HaplotypeCaller (GVCF mode)
    # ===========================================================
    # Calls variants per sample against the supercontig reference.
    # -ERC GVCF records confidence at every site (variant and
    # reference), required for joint genotyping in script 05.
    # Runtime is much shorter than the genome pipeline (~minutes
    # vs ~17h) because the reference is only ~1.9 Mb.

    echo "[$sample] Running HaplotypeCaller..."

    gatk --java-options "-Xmx16g" HaplotypeCaller \
        -R "$ref" \
        -I "$out/md/${sample}.sc.md.bam" \
        -O "$out/gvcf/${sample}.sc.g.vcf.gz" \
        -ERC GVCF \
        --native-pair-hmm-threads "$threads_hc" \
        2> "$out/logs/${sample}.hc.log"

    echo "[$sample] Done → $out/gvcf/${sample}.sc.g.vcf.gz"
}

export -f process_sample

# ============================================================
# --- Run in parallel -----------------------------------------
# ============================================================

ALL_R1=("$IN"/Gymno_*_1.trim.fq.gz)

echo "============================================================"
echo "STEP 04: Per-sample alignment and variant calling"
echo "Samples: ${#ALL_R1[@]}"
echo "Running $N_JOBS samples in parallel, $THREADS_BWA BWA threads each."
echo "Reference: $REF"
echo "Output: $OUT"
echo "============================================================"

parallel -j "$N_JOBS" \
    process_sample {} "$IN" "$OUT" "$REF" "$THREADS_BWA" "$THREADS_HC" \
    ::: "${ALL_R1[@]}"

echo ""
echo "============================================================"
echo "04_variantcall.sh complete."
echo "Key outputs:"
echo "  Alignment stats : $OUT/stats/*.sc.flagstat.txt"
echo "  Dup metrics     : $OUT/metrics/*.sc.metrics.txt"
echo "  Per-sample gVCFs: $OUT/gvcf/*.sc.g.vcf.gz"
echo "Next step: run 05_joint_genotyping.sh"
echo "============================================================"
