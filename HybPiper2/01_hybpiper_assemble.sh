#!/bin/bash
# ============================================================
# 01_hybpiper_assemble.sh — Target recovery with HybPiper2
# ============================================================
# This script treats whole-genome sequencing reads as if they
# came from a hybridization capture experiment, and recovers
# sequences for target loci defined by a probe set.
#
# WHAT HYBPIPER DOES (per sample):
# ---------------------------------
# 1. Maps all trimmed reads against the target probe sequences
#    using BWA. Most reads from WGS will NOT map (they come from
#    non-target regions of the genome) — only those that hit a
#    target locus are kept. This is the "virtual capture" step.
# 2. For each target locus, assembles the mapped reads into
#    contigs using SPAdes.
# 3. Identifies the contig most likely to be orthologous to the
#    target (i.e., the real gene copy, not a paralog).
# 4. Extracts the target region from the best contig.
# 5. Optionally recovers flanking intron sequence (--run_intronerate)
#    to produce "supercontig" = exon + surrounding intron.
#    Intron sequence has more variation and is valuable for
#    population-level analyses within a species.
#
# PROBE SET USED:
# ---------------
# Fabaceae1005_3273probes_1005reg.fasta
#   - 3,273 probe sequences targeting 1,005 conserved nuclear
#     loci distributed across the legume family (Fabaceae).
#   - Designed by Crameri et al. (2022, Mol. Ecol. Resources)
#     as a broadly applicable set validated across all six
#     Fabaceae subfamilies, including Caesalpinioideae, which
#     contains Gymnocladus.
#   - Starting with this set to assess how many loci recover
#     usable sequence in KCT before designing a species-specific
#     probe set or supplementing with additional loci.
#
# OUTPUT (one subdirectory per sample inside $OUT):
#   {prefix}/                 → HybPiper working directory
#   {prefix}/{locus}/         → per-locus assembly files
#   {prefix}/{prefix}.fasta   → assembled exon sequences
#   seq_lengths.tsv           → recovery stats table
#   recovery_heatmap.png      → heatmap of locus recovery across samples
#   retrieved_sequences/      → one multi-FASTA per locus (analysis-ready)
#
# PARALLELISM:
# ------------
# GNU parallel runs N_JOBS samples at the same time.
# Each sample uses THREADS CPU threads for BWA mapping.
# Total CPU = N_JOBS × THREADS.
# NOTE: WGS data is large (~5–11 GB per sample). Mapping all
# reads against the probe set takes longer than true capture
# data (which has far fewer off-target reads). Plan for ~30–60
# min per sample depending on coverage depth and core count.
#
# RUNTIME ESTIMATE:
# ~30-60 min per sample; with 8 parallel jobs → ~1-2 hours total
#
# Requirements: conda environment "hybpiper" must be active
#   source /home/nlove/miniconda3/etc/profile.d/conda.sh
#   conda activate hybpiper
# ============================================================

# Exit immediately if any command fails, if an unset variable
# is used, or if any command in a pipe fails
set -euo pipefail

# ============================================================
# --- Configuration -------------------------------------------
# ============================================================

# Trimmed FASTQ files (output of 03_fastp.sh)
IN="/home/nlove/kct_genomics/output_files/fastp/trimmed"

# Target FASTA file — 1,005 reference sequences from the Cajanus cajan
# genome, reformatted for HybPiper2 by 00_prep_targets.sh.
# Headers are in the required >GENEID-taxon format (e.g. >7271-Cajanus_cajan).
# See 00_prep_targets.sh for full explanation of what these sequences are.
TARGET="/home/nlove/kct_genomics/scripts/Fabaceae_probes_targets/Fabaceae_iter2_1005reg_hybpiper.fasta"

# Output directory — HybPiper creates one subfolder per sample here
OUT="/home/nlove/kct_genomics/output_files/hybpiper"

# --- Thread and parallelism settings -------------------------
N_JOBS=8       # Number of samples to process simultaneously
               # Adjust based on available cores: N_JOBS × THREADS = total cores used
THREADS=2      # CPU threads per sample (used by BWA for read mapping)
               # 8 jobs × 2 threads = 16 cores total (test run: 5 × 2 = 10 cores)

# --- Test mode -----------------------------------------------
# Set TEST_MODE=true to run only the first 5 samples as a quick
# check before committing to the full 61-sample run.
# Confirms: HybPiper is working, target file is compatible,
# recovery looks reasonable, and runtime is as expected.
# Set to false (or remove) once you're happy with the test results.
TEST_MODE=true

# ============================================================
# --- Setup ---------------------------------------------------
# ============================================================

mkdir -p "$OUT"

# HybPiper creates output directories in the CURRENT working directory.
# We must cd to $OUT before running assemble so all sample folders
# are created there (not wherever this script was launched from).
cd "$OUT"

# ============================================================
# --- STEP 1: Assemble (map reads → assemble per locus) -------
# ============================================================
# This is the main HybPiper step and the most time-consuming.
# Each sample is processed independently, so we parallelize
# across samples using GNU parallel.

# Build the list of R1 files to process
ALL_R1=("$IN"/*_1.trim.fq.gz)

if [[ "$TEST_MODE" == true ]]; then
    # Slice to first 5 samples for the test run
    R1_FILES=("${ALL_R1[@]:0:5}")
    echo "TEST MODE: running first 5 samples only."
else
    R1_FILES=("${ALL_R1[@]}")
fi

echo "============================================================"
echo "STEP 1: HybPiper assemble"
echo "Samples to process: ${#R1_FILES[@]} (of ${#ALL_R1[@]} total)"
echo "Running $N_JOBS samples in parallel, $THREADS threads each."
echo "Target file: $TARGET"
echo "Output directory: $OUT"
echo "============================================================"

# Define the function that processes one sample.
# GNU parallel will call this once per R1 file.
assemble_sample() {
    local r1="$1"
    local in="$2"
    local target="$3"
    local threads="$4"
    local out="$5"

    # Construct the R2 filename (replace _1 with _2)
    local r2="${r1/_1.trim.fq.gz/_2.trim.fq.gz}"

    # Safety check: skip if R2 is missing
    if [[ ! -f "$r2" ]]; then
        echo "WARNING: missing R2 for $r1 — skipping"
        return 0
    fi

    # Extract full base name (without path or _1.trim.fq.gz suffix)
    # e.g. Gymno_2013_109_CKDL250009712-1A_22V7CCLT4_L2
    local base
    base=$(basename "$r1" _1.trim.fq.gz)

    # --- Skip Undetermined demultiplexing artifacts ------------
    if [[ "$base" == Undetermined* ]]; then
        echo "Skipping $base — Undetermined reads, not a real sample"
        return 0
    fi

    # Clean sample name: keep only Gymno_YYYY_NNN (first 3 fields)
    # This becomes the output directory name and sample ID in all
    # downstream outputs (stats tables, alignments, tree tips).
    # e.g. Gymno_2013_109_CKDL250009712-1A_22V7CCLT4_L2 → Gymno_2013_109
    local prefix
    prefix=$(echo "$base" | cut -d'_' -f1-3)

    # --- Resume/skip check ------------------------------------
    # If this sample's output directory already exists and contains
    # assembled sequences, skip it. Safe to re-run after interruption.
    if [[ -f "$out/${prefix}/${prefix}.fasta" ]]; then
        echo "Skipping $prefix — assembly already exists"
        return 0
    fi

    echo "Assembling: $prefix"

    # --- Run HybPiper assemble --------------------------------
    # Note: must be run from $OUT (cd above) so the {prefix}/ output
    # directory is created in the right place.
    #
    # -t_dna $TARGET
    #   Target sequences as nucleotide FASTA. HybPiper maps reads
    #   against these, then assembles only the reads that mapped.
    #
    # -r $r1 $r2
    #   Paired trimmed FASTQ files for this sample. No --unpaired
    #   flag needed here — fastp produces only paired output.
    #
    # --prefix $prefix
    #   Output directory name and sample label (Gymno_YYYY_NNN).
    #
    # --bwa
    #   Use BWA for mapping (required when target is nucleotide, -t_dna).
    #   Alternative is DIAMOND (for protein targets, -t_aa) — not used here.
    #
    # --run_intronerate
    #   After assembling exon sequences, also recover flanking intron
    #   sequence ("supercontigs" = exon + intron). Introns have higher
    #   variation within a species and are valuable for population-level
    #   analyses. Note: HybPiper ≥ 2.1.6 runs this by default; explicit
    #   flag included here for clarity and compatibility with older versions.
    #
    # --cpu $threads
    #   Number of threads for BWA mapping within this sample.
    hybpiper assemble \
        -t_dna "$target" \
        -r "$r1" "$r2" \
        --prefix "$prefix" \
        --bwa \
        --run_intronerate \
        --cpu "$threads" \
        2>&1 | tee "$out/${prefix}.assemble.log"

    echo "Done: $prefix"
}

# Export the function so GNU parallel can access it in subprocesses
export -f assemble_sample

# Run assemble_sample in parallel across all R1 files
parallel -j "$N_JOBS" \
    assemble_sample {} "$IN" "$TARGET" "$THREADS" "$OUT" \
    ::: "${R1_FILES[@]}"

echo ""
echo "STEP 1 complete. Sample directories in: $OUT"

# ============================================================
# --- STEP 2: Generate a sample name file ---------------------
# ============================================================
# hybpiper stats and retrieve_sequences need a plain text file
# listing one sample prefix (Gymno_YYYY_NNN) per line.

echo "============================================================"
echo "STEP 2: Generating sample name file"
echo "============================================================"

NAMEFILE="$OUT/namefile.txt"

# Build the name file from the assembled sample directories.
# Each line = the prefix used in STEP 1 (Gymno_YYYY_NNN).
# We find these by looking for the characteristic .fasta output file.
for fasta in "$OUT"/Gymno_*/*.fasta; do
    basename "$(dirname "$fasta")"
done > "$NAMEFILE"

echo "$(wc -l < "$NAMEFILE") samples written to $NAMEFILE"

# ============================================================
# --- STEP 3: Summary statistics ------------------------------
# ============================================================
# Produces seq_lengths.tsv: a table showing, for each sample × locus,
# how much of the target sequence was recovered (as % of target length).
# Values close to 1.0 = full recovery; 0 = nothing recovered.
# This table is also the input for the recovery heatmap (STEP 4).

echo "============================================================"
echo "STEP 3: hybpiper stats"
echo "============================================================"

hybpiper stats \
    -t_dna "$TARGET" \
    gene \
    "$NAMEFILE"

echo "Stats complete. Output: $OUT/seq_lengths.tsv"

# ============================================================
# --- STEP 4: Recovery heatmap --------------------------------
# ============================================================
# Generates a PNG heatmap: rows = samples, columns = loci.
# Color intensity = fraction of target sequence recovered.
# This is your first diagnostic — it tells you how many loci
# are usable across your samples for downstream analyses.

echo "============================================================"
echo "STEP 4: hybpiper recovery_heatmap"
echo "============================================================"

hybpiper recovery_heatmap seq_lengths.tsv

echo "Heatmap written to: $OUT/recovery_heatmap.png"

# ============================================================
# --- STEP 5: Retrieve sequences ------------------------------
# ============================================================
# Collects the assembled target sequences for all samples into
# one multi-FASTA file per locus — the format needed for
# alignment and population-level analyses.
#
# "dna" retrieves exon sequences only.
# To retrieve exon + intron (supercontig), replace "dna" with
# "supercontig" — but only if --run_intronerate ran successfully.

echo "============================================================"
echo "STEP 5: hybpiper retrieve_sequences"
echo "============================================================"

hybpiper retrieve_sequences dna \
    -t_dna "$TARGET" \
    --sample_names "$NAMEFILE"

echo "Sequences retrieved."
echo ""
echo "============================================================"
echo "HybPiper2 pipeline complete."
echo "Key outputs:"
echo "  Recovery heatmap : $OUT/recovery_heatmap.png"
echo "  Stats table      : $OUT/seq_lengths.tsv"
echo "  Retrieved seqs   : $OUT/*.FNA  (one file per locus)"
echo "============================================================"
