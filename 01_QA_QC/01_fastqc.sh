#!/bin/bash
# Requirements: conda environment "qc" must be active
#   (contains fastqc, fastp, multiqc)
#   source /home/nlove/miniconda3/etc/profile.d/conda.sh
#   conda activate qc

BASE="/home/nlove/kct_genomics"
IN="$BASE/input_files/NovaSeq/01.RawData/all_fastq"
OUT="$BASE/output_files/fastqc_results"

THREADS=20

mkdir -p "$OUT"

fastqc -t "$THREADS" -o "$OUT" "$IN"/*.fq.gz
