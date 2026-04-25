#!/bin/bash
# ============================================================
# 02_paralog_retriever.sh — Identify and summarize paralogous loci
# ============================================================
# Runs hybpiper paralog_retriever across all 61 samples to produce
# a cross-sample summary of which loci show paralog signal.
#
# WHAT ARE PARALOGS IN THIS CONTEXT?
# -----------------------------------
# After HybPiper assembles contigs for each locus, it flags two
# types of potential paralogy (both recorded in hybpiper_stats.tsv):
#
#   ParalogWarningsLong:
#     SPAdes assembled multiple contigs for this locus, and at
#     least one extra contig is ≥75% the length of the reference.
#     This suggests a second gene copy was captured alongside the
#     intended target.
#
#   ParalogWarningsDepth:
#     Read depth across the locus shows two distinct coverage
#     peaks — the signature of two similar sequences being mapped
#     together. This indicates the locus may be duplicated in the
#     Gymnocladus genome.
#
# Across your 61 samples, ~50–71 loci per sample trigger the Long
# warning and ~86–125 trigger the Depth warning (~150–190 unique
# loci per sample with some paralog signal, ~15–19% of 1005).
#
# WHAT DOES paralog_retriever DO?
# ---------------------------------
# It reads the per-sample paralog FASTA files written by
# hybpiper assemble and produces a cross-sample summary:
#
#   paralogs_all/
#     One FASTA per flagged locus, containing all paralog
#     sequences recovered across all samples. Used if you want
#     to inspect or align the duplicated sequences directly.
#
#   paralog_heatmap.png
#     Heatmap (loci × samples) showing which loci are flagged
#     for paralogy in which samples. Loci flagged in many or
#     all samples are likely true multi-copy genes in Gymnocladus
#     and should be excluded from population genomics.
#
#   putative_paralog_report.tsv (or similar)
#     Table summarising how many samples flag each locus.
#     Use this to apply a per-locus exclusion cutoff:
#       - Conservative: exclude any locus flagged in ≥1 sample
#         → expect ~800+ clean loci
#       - Lenient: exclude only loci flagged in the majority
#         of samples
#
# NEXT STEPS AFTER THIS SCRIPT
# -----------------------------
# 1. Inspect paralog_heatmap.png — look for loci with broad
#    signal across many samples vs. sample-specific flags.
# 2. Decide on exclusion cutoff and generate a clean locus list.
# 3. Proceed to scripts 03+ for alignment and variant calling
#    on the clean locus set.
#
# Requirements: conda environment "hybpiper" must be active
#   source /home/nlove/miniconda3/etc/profile.d/conda.sh
#   conda activate hybpiper
#
# Must be run from $OUT (the hybpiper assemble output directory)
# so paralog_retriever can find per-sample subdirectories.
# ============================================================

set -euo pipefail

# ============================================================
# --- Configuration -------------------------------------------
# ============================================================

# Target FASTA file (same as used in 01_hybpiper_assemble.sh)
TARGET="/home/nlove/kct_genomics/scripts/Fabaceae_probes_targets/Fabaceae_iter2_1005reg_hybpiper.fasta"

# HybPiper assemble output directory — contains one subdirectory
# per sample. paralog_retriever must be run from here.
OUT="/home/nlove/kct_genomics/output_files/hybpiper"

# Sample name file generated in 01_hybpiper_assemble.sh STEP 2
NAMEFILE="$OUT/namefile.txt"

# ============================================================
# --- Run paralog_retriever -----------------------------------
# ============================================================

cd "$OUT"

echo "============================================================"
echo "hybpiper paralog_retriever"
echo "Target file : $TARGET"
echo "Namefile    : $NAMEFILE"
echo "Output dir  : $OUT"
echo "============================================================"

hybpiper paralog_retriever \
    --targetfile_dna "$TARGET" \
    --namelist "$NAMEFILE" \
    --heatmap_dpi 300

echo ""
echo "============================================================"
echo "paralog_retriever complete."
echo "Key outputs:"
echo "  Paralog heatmap  : $OUT/paralog_heatmap.png"
echo "  Paralog sequences: $OUT/paralogs_all/"
echo "  Summary table    : $OUT/putative_paralog_report.tsv"
echo "============================================================"
