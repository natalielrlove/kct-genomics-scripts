#!/bin/bash

BASE="/home/nlove/kct_genomics"
IN="$BASE/input_files/NovaSeq/01.RawData/all_fastq"
OUT="$BASE/output_files/fastp"

mkdir -p "$OUT/trimmed" "$OUT/html" "$OUT/json"

for r1 in "$IN"/*_1.fq.gz; do
  r2="${r1%_1.fq.gz}_2.fq.gz"

  if [ ! -f "$r2" ]; then
    echo "Skipping: missing R2 for $r1"
    continue
  fi

  base=$(basename "$r1" _1.fq.gz)

  fastp \
    -i "$r1" -I "$r2" \
    -o "$OUT/trimmed/${base}_1.trim.fq.gz" \
    -O "$OUT/trimmed/${base}_2.trim.fq.gz" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --thread 6 \
    --html "$OUT/html/${base}.fastp.html" \
    --json "$OUT/json/${base}.fastp.json"
done
