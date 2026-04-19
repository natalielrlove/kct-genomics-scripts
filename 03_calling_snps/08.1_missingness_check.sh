#!/bin/bash
# ============================================================
# 08.1_missingness_check.sh — Evaluate missingness thresholds
# ============================================================
# Checks how many PASS SNPs are retained at different levels
# of per-site missingness (F_MISSING). Helps choose a threshold
# for population-level filtering in script 10.
#
# F_MISSING is the fraction of samples with no genotype call
# at a site (./.). A site with F_MISSING = 0.3 means 30% of
# samples (i.e. ~18 of 61) have missing data there.
#
# Run the chr01-only version first (~minutes) to get a fast
# preview, then run the full VCF (~20-30 min) to confirm.
#
# Requirements: conda environment "gatk" must be active
#   (contains bcftools)
# ============================================================

set -euo pipefail

# Activate conda environment before running this script:
#   source /home/nlove/miniconda3/etc/profile.d/conda.sh
#   conda activate gatk

VCF=~/kct_genomics/output_files/gatk_filter_hap1/gymno_hap1.snps.filtered.vcf.gz

# --- chr01 preview (fast) ------------------------------------
echo "=== chr01 only (preview) ==="
total_chr01=$(bcftools stats -r chr01 "$VCF" | awk '/^SN.*number of records/ {print $NF}')
for miss in 0.5 0.3 0.1; do
  kept=$(bcftools view -r chr01 -i "F_MISSING <= ${miss}" "$VCF" | bcftools view -H | wc -l)
  pct=$(echo "scale=1; ${kept} * 100 / ${total_chr01}" | bc)
  echo "F_MISSING <= ${miss}: ${kept} / ${total_chr01} sites kept (${pct}%)"
done

# --- Full VCF (slow, ~20-30 min) -----------------------------
echo ""
echo "=== Full VCF ==="
total=$(bcftools stats "$VCF" | awk '/^SN.*number of records/ {print $NF}')
for miss in 0.5 0.3 0.1; do
  kept=$(bcftools view -i "F_MISSING <= ${miss}" "$VCF" | bcftools view -H | wc -l)
  pct=$(echo "scale=1; ${kept} * 100 / ${total}" | bc)
  echo "F_MISSING <= ${miss}: ${kept} / ${total} sites kept (${pct}%)"
done
