#!/bin/bash
# ============================================================
# 07_genomicsdb_import.sh — Combine per-sample gVCFs and joint genotype
# ============================================================
# This script takes the 61 per-sample gVCF files produced by
# 06_haplotypecaller.sh and combines them into a single raw
# multi-sample VCF through two sequential steps:
#
#   Step 1 — Generate a sample map (input list of samples for GATK; look up table)
#   Step 2 — Generate an interval list from the reference index
#   Step 3 — GenomicsDBImport: load all 61 gVCFs into a database
#   Step 4 — GenotypeGVCFs: query the database to produce the final VCF
#
# WHY JOINT GENOTYPING?
# ---------------------
# Running HaplotypeCaller per sample (script 06) was intentional —
# it is accurate and parallelizable. But calling each sample alone
# means rare variants can look like sequencing errors (only one
# sample supports them). Joint genotyping fixes this by looking at
# all 61 samples together: if a rare variant is consistently seen
# in even one sample, and all other samples clearly show the
# reference at that site, GATK can confidently distinguish
# "real rare variant" from "error."
#
# WHY GenomicsDBImport instead of CombineGVCFs?
# ----------------------------------------------
# Both tools combine per-sample gVCFs before joint genotyping.
# CombineGVCFs merges them into one big gVCF file.
# GenomicsDBImport loads them into a GenomicsDB data store —
# a column-oriented database optimized for querying genomic data.
#
# GenomicsDBImport is the GATK-recommended approach because:
#   - Faster for large cohorts (reads all samples simultaneously
#     rather than sequentially merging files)
#   - Lower memory footprint at scale
#   - The resulting database is queried directly by GenotypeGVCFs
#     using the "gendb://" prefix (see Step 4)
#
# Think of it like this:
#   CombineGVCFs = stacking 61 spreadsheets into one giant spreadsheet
#   GenomicsDBImport = loading 61 spreadsheets into a database with an index
#
# THREADING NOTE
# --------------
# This step is mostly single-threaded. The only meaningful parallelism
# option is --reader-threads (reads input gVCFs in parallel before
# the single-threaded merge). We use 4 reader threads here.
# Total expected runtime: ~40 hours (based on collaborator benchmark).
#
# OUTPUT
# ------
# GenomicsDB workspace (Step 3):
#   output_files/gatk_joint_hap1/genomicsdb/
#   (a directory, not a file — this is the database)
#
# Raw joint VCF (Step 4):
#   output_files/gatk_joint_hap1/gymno_hap1.raw.vcf.gz
#   (unfiltered; all 61 samples as columns, one row per variant site)
#
# RESUME BEHAVIOR
# ---------------
# - GenomicsDBImport: will exit with a helpful message if the workspace
#   directory already exists (GATK cannot overwrite it automatically).
#   To re-run from scratch: rm -rf output_files/gatk_joint_hap1/genomicsdb/
# - GenotypeGVCFs: skipped if the output VCF already exists.
#
# Requirements: conda environment "gatk" must be active
#  see below
# ============================================================

# Exit immediately if any command fails, if an unset variable
# is used, or if any command in a pipe fails
set -euo pipefail

# --- Activate conda environment ------------------------------
# The "gatk" environment contains GATK4.
# Activate manually before running this script:
#   source /home/nlove/miniconda3/etc/profile.d/conda.sh
#   conda activate gatk


# --- Directory paths -----------------------------------------
# Reference genome (haplotype 1) — same reference used throughout pipeline
REF="/home/nlove/kct_genomics/input_files/reference/Gymnocladus_dioicus_M_hap1.fa"

# Input: per-sample gVCF files produced by 06_haplotypecaller.sh
GVCF_DIR="/home/nlove/kct_genomics/output_files/gatk_gvcf_hap1/gvcf"

# Output: all joint genotyping outputs go here
OUT="/home/nlove/kct_genomics/output_files/gatk_joint_hap1"

# --- Output file names ---------------------------------------
WORKSPACE="$OUT/genomicsdb"         # GenomicsDB data store (directory, not a file)
SAMPLE_MAP="$OUT/tmp/sample_map.tsv"    # Tab-separated list of sample name → gVCF path
INTERVALS="$OUT/tmp/all_scaffolds.list" # List of all scaffold names from the reference
RAW_VCF="$OUT/gymno_hap1.raw.vcf.gz"   # Final output: raw joint-genotyped VCF
LOG="$OUT/logs/07_genomicsdb_import.log"

# --- Create output directories if they don't already exist ---
# tmp/  = sample map, interval list, and GATK temporary files
# logs/ = log file for this script
mkdir -p "$OUT"/{logs,tmp}

# ============================================================
# STEP 1 — Generate sample map
# ============================================================
# GenomicsDBImport requires a sample map: a tab-separated file
# with one row per sample in the format:
#   sample_name <TAB> /path/to/sample.g.vcf.gz
#
# The sample name is extracted from the gVCF filename by stripping
# the directory path and the .g.vcf.gz extension.
# e.g.: Gymno_2013_109_CKDL250009712-1A_22V7CCLT4_L2.g.vcf.gz
#     → Gymno_2013_109_CKDL250009712-1A_22V7CCLT4_L2  (sample name)
#
# This is analogous to the sample sheet you might create for other tools.

echo "=== Step 1: Generating sample map ===" | tee "$LOG"
echo "Started: $(date)" | tee -a "$LOG"

for gvcf in "$GVCF_DIR"/*.g.vcf.gz; do
    sample=$(basename "$gvcf" .g.vcf.gz)
    printf "%s\t%s\n" "$sample" "$gvcf"
done > "$SAMPLE_MAP"

echo "Sample map written to: $SAMPLE_MAP" | tee -a "$LOG"
echo "Samples in map: $(wc -l < "$SAMPLE_MAP")" | tee -a "$LOG"

# ============================================================
# STEP 2 — Generate interval list
# ============================================================
# GenomicsDBImport requires -L (intervals) — it will not run without it.
# Intervals tell the tool which regions of the genome to include.
# We want ALL regions, so we include all scaffolds.
#
# The reference FASTA index (.fai file) has one row per scaffold:
#   scaffold_name <TAB> length <TAB> offset <TAB> ...
# We extract column 1 (scaffold names) to build a simple list.
# With 986 scaffolds in this assembly, this list will have 986 lines.

echo "" | tee -a "$LOG"
echo "=== Step 2: Generating interval list ===" | tee -a "$LOG"

awk '{print $1}' "${REF}.fai" > "$INTERVALS"

echo "Interval list written to: $INTERVALS" | tee -a "$LOG"
echo "Scaffolds in interval list: $(wc -l < "$INTERVALS")" | tee -a "$LOG"

# ============================================================
# STEP 3 — GenomicsDBImport
# ============================================================
# Loads all 61 per-sample gVCFs into a GenomicsDB data store.
# The result is a directory ("workspace") containing a compressed,
# indexed database of all variant and reference-confidence records
# across all samples. It is NOT a VCF file — it is a database
# that GenotypeGVCFs (Step 4) will query directly.
#
# This step is mostly single-threaded and memory-intensive.
# Expected runtime: several hours (exact depends on server load).

echo "" | tee -a "$LOG"
echo "=== Step 3: GenomicsDBImport ===" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"

# Resume check: GenomicsDBImport will fail if the workspace directory
# already exists (returns TRUE). Exit early with a helpful message rather than
# letting GATK fail cryptically.
if [[ -d "$WORKSPACE" ]]; then
    echo "ERROR: GenomicsDB workspace already exists at:" | tee -a "$LOG"
    echo "  $WORKSPACE" | tee -a "$LOG"
    echo "GenomicsDBImport cannot overwrite an existing workspace." | tee -a "$LOG"
    echo "To re-run from scratch, delete it first:" | tee -a "$LOG"
    echo "  rm -rf $WORKSPACE" | tee -a "$LOG"
    exit 1
fi

gatk --java-options "-Xmx64g -Djava.io.tmpdir=$OUT/tmp" GenomicsDBImport \
    --sample-name-map "$SAMPLE_MAP" \
    --genomicsdb-workspace-path "$WORKSPACE" \
    -L "$INTERVALS" \
    --merge-input-intervals \
    --reader-threads 4 \
    --tmp-dir "$OUT/tmp" \
    2>&1 | tee -a "$LOG"
#
# --java-options "-Xmx64g"
#   Allow Java to use up to 64 GB RAM. GenomicsDBImport holds data
#   from all 61 samples in memory during import — more than
#   HaplotypeCaller needed per sample. Curie has 481 GB available
#   so 64 GB is conservative but comfortable.
#
# --sample-name-map
#   The tab-separated file from Step 1, mapping each sample name
#   to its gVCF path. GATK reads this instead of requiring one
#   -V flag per sample on the command line (which would be unwieldy
#   for 61 samples).
#
# --genomicsdb-workspace-path
#   Where to create the GenomicsDB database directory.
#   This directory must NOT already exist — GATK creates it fresh.
#
# -L (passed as the interval list file)
#   All 986 scaffolds. Without this, GATK refuses to run.
#
# --merge-input-intervals
#   With 986 scaffolds, GATK would by default create a separate
#   database partition per scaffold. This flag merges them all into
#   one database, which is simpler and avoids per-scaffold overhead.
#   Recommended for fragmented assemblies like this one.
#
# --reader-threads 4
#   The only real parallelism available here: 4 threads read the
#   input gVCF files in parallel before the single-threaded merge.
#   More than 4-6 provides diminishing returns.
#
# --tmp-dir
#   GenomicsDBImport writes large temporary files during import.
#   Point this to a project directory with plenty of space —
#   the default /tmp can fill up with 61 samples.

echo "" | tee -a "$LOG"
echo "GenomicsDBImport complete: $(date)" | tee -a "$LOG"
echo "Database written to: $WORKSPACE" | tee -a "$LOG"

# ============================================================
# STEP 4 — GenotypeGVCFs
# ============================================================
# Queries the GenomicsDB database and performs joint genotyping
# across all 61 samples. For every site in the genome where at
# least one sample shows evidence of a variant, GATK computes
# genotype likelihoods for all samples simultaneously and outputs
# the result as a standard multi-sample VCF.
#
# The output is "raw" — variants are present but not yet filtered.
# Filtering (VQSR or hard filters) is the next pipeline step.
#
# This step is also single-threaded and memory-intensive.
# Expected runtime: several hours (runs after GenomicsDBImport completes).

echo "" | tee -a "$LOG"
echo "=== Step 4: GenotypeGVCFs ===" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"

# Resume check: skip genotyping if the output VCF already exists
if [[ -f "$RAW_VCF" ]]; then
    echo "Raw VCF already exists: $RAW_VCF — skipping GenotypeGVCFs." | tee -a "$LOG"
    exit 0
fi

gatk --java-options "-Xmx64g -Djava.io.tmpdir=$OUT/tmp" GenotypeGVCFs \
    -R "$REF" \
    -V "gendb://$WORKSPACE" \
    -O "$RAW_VCF" \
    --tmp-dir "$OUT/tmp" \
    2>&1 | tee -a "$LOG"
#
# --java-options "-Xmx64g"
#   GenotypeGVCFs also needs substantial memory — it reads the
#   entire database and computes likelihoods across all 61 samples.
#
# -R : reference genome FASTA
#
# -V "gendb://$WORKSPACE"
#   The "gendb://" prefix tells GATK to read from a GenomicsDB
#   data store rather than a VCF file. This is how GenotypeGVCFs
#   connects to the database created in Step 3.
#   Without the prefix it would expect a plain VCF path.
#
# -O : output raw joint VCF (compressed .vcf.gz; .tbi index created automatically)
#   This is the final product of joint genotyping: one VCF with
#   all 61 samples as columns and one row per variant site.
#
# --tmp-dir : same as above, avoid filling /tmp

echo "" | tee -a "$LOG"
echo "==================================================" | tee -a "$LOG"
echo "=== All steps complete ===" | tee -a "$LOG"
echo "Finished: $(date)" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "GenomicsDB workspace : $WORKSPACE" | tee -a "$LOG"
echo "Raw joint VCF        : $RAW_VCF" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "Next step: variant filtering (script 08)" | tee -a "$LOG"
