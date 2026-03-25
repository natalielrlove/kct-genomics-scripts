#!/bin/bash
# Requirements: conda environment "qc" must be active
#   (contains fastqc, fastp, multiqc)
#   source /home/nlove/miniconda3/etc/profile.d/conda.sh
#   conda activate qc

BASE="/home/nlove/kct_genomics"
IN="$BASE/output_files/fastqc_results"
OUT="$BASE/output_files/multiqc"

mkdir -p "$OUT"

multiqc "$IN" -o "$OUT" --force
