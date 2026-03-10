#!/bin/bash

BASE="/home/nlove/kct_genomics"
IN="$BASE/output_files/fastqc_results"
OUT="$BASE/output_files/multiqc"

mkdir -p "$OUT"

multiqc "$IN" -o "$OUT" --force
