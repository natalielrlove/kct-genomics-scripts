#!/bin/bash
# ============================================================
# 03_fastp.sh — Adapter trimming and quality filtering
# ============================================================
# This script runs fastp on all raw paired-end FASTQ files.
# fastp does two things:
#   1. Removes adapter sequences from the ends of reads
#   2. Filters out low-quality reads
# For each sample it produces:
#   - Two trimmed FASTQ files (R1 and R2), ready for alignment
#   - An HTML report (visual summary of QC metrics)
#   - A JSON report (machine-readable, used by MultiQC)
# ============================================================

# --- Directory paths -----------------------------------------
# BASE: root directory for all project files on the server
BASE="/home/nlove/kct_genomics"

# IN: folder containing all raw FASTQ files from the sequencer
# Files are named: SampleName_1.fq.gz (R1) and SampleName_2.fq.gz (R2)
IN="$BASE/input_files/NovaSeq/01.RawData/all_fastq"

# OUT: folder where all fastp outputs will be written
OUT="$BASE/output_files/fastp"

# --- Create output directories if they don't already exist ---
# -p means: create parent directories as needed, no error if they exist
# trimmed/ = cleaned FASTQ files
# html/    = per-sample visual QC reports
# json/    = per-sample machine-readable reports (used by MultiQC)
mkdir -p "$OUT/trimmed" "$OUT/html" "$OUT/json"

# --- Loop over every R1 (forward read) file ------------------
# The wildcard *_1.fq.gz matches one file per sample (the forward read).
# We process both R1 and R2 together for each sample.
for r1 in "$IN"/*_1.fq.gz; do

  # Construct the R2 filename by replacing _1.fq.gz with _2.fq.gz
  # e.g. Gymno_2013_109_..._1.fq.gz → Gymno_2013_109_..._2.fq.gz
  r2="${r1%_1.fq.gz}_2.fq.gz"

  # Safety check: skip this sample if the R2 file is missing
  # (would indicate a problem with the sequencing delivery)
  if [ ! -f "$r2" ]; then
    echo "Skipping: missing R2 for $r1"
    continue
  fi

  # Extract just the sample name (no directory path, no _1.fq.gz suffix)
  # e.g. /path/to/Gymno_2013_109_CKDL..._1.fq.gz → Gymno_2013_109_CKDL...
  # This is used to name all output files consistently
  base=$(basename "$r1" _1.fq.gz)

  # --- Run fastp -----------------------------------------------
  fastp \
    -i "$r1" -I "$r2" \                                        # Input: R1 and R2 raw FASTQ files
    -o "$OUT/trimmed/${base}_1.trim.fq.gz" \                   # Output: trimmed R1
    -O "$OUT/trimmed/${base}_2.trim.fq.gz" \                   # Output: trimmed R2
    --detect_adapter_for_pe \                                  # Auto-detect adapter sequences for paired-end data
    --qualified_quality_phred 20 \                             # Discard bases with Phred quality score < 20 (~99% accuracy)
    --length_required 50 \                                     # Discard reads shorter than 50 bp after trimming
    --thread 20 \                                              # Use 20 CPU threads to run faster
    --html "$OUT/html/${base}.fastp.html" \                    # Write per-sample HTML report here
    --json "$OUT/json/${base}.fastp.json"                      # Write per-sample JSON report here (used by MultiQC)

done
