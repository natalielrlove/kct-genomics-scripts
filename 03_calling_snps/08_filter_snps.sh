#!/bin/bash
# ============================================================
# 08_filter_snps.sh — Apply GATK hard filters to the raw joint VCF
# ============================================================
# This script takes the raw multi-sample VCF produced by
# 07_genomicsdb_import.sh and produces a filtered SNP-only VCF
# through three sequential steps:
#
#   Step 1 — SelectVariants: extract SNPs only (discard indels)
#   Step 2 — VariantFiltration: tag SNPs that fail quality thresholds
#   Step 3 — SelectVariants: remove tagged (FILTER != PASS) sites
#
# WHY HARD FILTERS?
# -----------------
# GATK's preferred filtering method is VQSR (Variant Quality Score
# Recalibration), which trains a statistical model to distinguish
# true variants from artifacts. However, VQSR requires a large set
# of known variants (e.g., dbSNP, HapMap) as a truth set — these
# resources only exist for humans and a handful of model organisms.
# For Gymnocladus dioicus, no such truth set exists.
#
# Hard filters are the GATK-recommended alternative for non-model
# organisms. Instead of a trained model, they apply fixed thresholds
# to annotations that GATK computes during variant calling. Sites
# that fall outside the expected range for real variants are flagged
# as likely technical artifacts and removed.
#
# WHY THREE STEPS?
# ----------------
# Step 1 (SelectVariants) separates SNPs from indels before
# filtering because GATK recommends different hard filter thresholds
# for each variant type. This script handles SNPs only.
#
# Step 2 (VariantFiltration) does NOT remove sites — it adds a
# FILTER field to each site recording which thresholds it failed.
# Sites that pass all filters get FILTER=PASS.
# Sites that fail get FILTER=<reason> (e.g., FILTER=FS60).
# This tagging-before-removal approach lets you inspect which
# filters are removing the most sites before committing.
#
# Step 3 (SelectVariants --exclude-filtered) removes all tagged
# sites, keeping only FILTER=PASS sites for downstream analysis.
#
# HARD FILTER THRESHOLDS (GATK recommended for SNPs)
# ---------------------------------------------------
# QD   < 2.0   QualByDepth: variant quality divided by depth.
#              Low QD means the variant call is weak relative to
#              how many reads support it — likely a low-confidence call.
#
# QUAL < 30.0  Phred-scaled confidence that the variant exists at this site.
#              QUAL 30 = 1-in-1000 error probability. Low QUAL means GATK
#              has weak overall confidence the variant is real, independent
#              of read depth (unlike QD).
#
# FS  > 60.0   FisherStrand: strand bias measured by Fisher's exact test.
#              High FS means the variant is seen predominantly on one
#              strand — a common sequencing artifact.
#
# SOR > 3.0    StrandOddsRatio: a more robust strand bias metric.
#              Complements FS; high SOR also indicates strand bias.
#
# MQ  < 40.0   RMSMappingQuality: average mapping quality of reads
#              covering the site. Low MQ means reads are mapping
#              poorly, suggesting a repetitive or ambiguous region.
#
# MQRankSum < -12.5   MappingQualityRankSum: tests whether the
#              mapping quality of reads supporting the alt allele
#              differs from reads supporting the ref allele.
#              A large negative value suggests the alt allele is
#              supported only by poorly-mapping reads (artifact signal).
#
# ReadPosRankSum < -8.0   ReadPosRankSum: tests whether the alt
#              allele tends to appear near the ends of reads.
#              Variants at read ends are more likely to be errors.
#
# OUTPUT
# ------
# SNPs extracted from raw VCF (before filtering):
#   output_files/gatk_filter_hap1/gymno_hap1.snps.vcf.gz
#
# SNPs with filter tags applied (for inspection):
#   output_files/gatk_filter_hap1/gymno_hap1.snps.tagged.vcf.gz
#
# Final filtered SNPs — PASS sites only (input for script 09):
#   output_files/gatk_filter_hap1/gymno_hap1.snps.filtered.vcf.gz
#
# RESUME BEHAVIOR
# ---------------
# Each output file is checked before running its step.
# Safe to re-run after interruption — completed steps are skipped.
#
# Requirements: conda environment "gatk" must be active
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

# Input: raw joint VCF produced by 07_genomicsdb_import.sh
RAW_VCF="/home/nlove/kct_genomics/output_files/gatk_joint_hap1/gymno_hap1.raw.vcf.gz"

# Output: all filtering outputs go here
OUT="/home/nlove/kct_genomics/output_files/gatk_filter_hap1"

# --- Output file names ---------------------------------------
SNP_VCF="$OUT/gymno_hap1.snps.vcf.gz"           # Step 1: SNPs only
TAGGED_VCF="$OUT/gymno_hap1.snps.tagged.vcf.gz" # Step 2: filter tags applied
FILTERED_VCF="$OUT/gymno_hap1.snps.filtered.vcf.gz" # Step 3: PASS sites only
LOG="$OUT/logs/08_filter_snps.log"

# --- Create output directories if they don't already exist ---
mkdir -p "$OUT"/{logs,tmp}

echo "=== 08_filter_snps.sh ===" | tee "$LOG"
echo "Started: $(date)" | tee -a "$LOG"
echo "Input VCF: $RAW_VCF" | tee -a "$LOG"

# ============================================================
# STEP 1 — SelectVariants: extract SNPs only
# ============================================================
# The raw VCF contains both SNPs and indels (insertions/deletions).
# GATK recommends filtering SNPs and indels separately because
# the hard filter thresholds differ between variant types.
# This step extracts only biallelic SNPs for filtering.
#
# Multiallelic sites (more than one alt allele) are excluded here
# with --restrict-alleles-to BIALLELIC. Most population genomics
# tools assume biallelic variants, and multiallelic SNPs require
# special handling that is out of scope for this pipeline.

echo "" | tee -a "$LOG"
echo "=== Step 1: Extracting SNPs ===" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"

if [[ -f "$SNP_VCF" ]]; then
    echo "SNP VCF already exists: $SNP_VCF — skipping Step 1." | tee -a "$LOG"
else
    gatk --java-options "-Xmx16g -Djava.io.tmpdir=$OUT/tmp" SelectVariants \
        -R "$REF" \
        -V "$RAW_VCF" \
        -O "$SNP_VCF" \
        --select-type-to-include SNP \
        --restrict-alleles-to BIALLELIC \
        2>&1 | tee -a "$LOG"
    #
    # --select-type-to-include SNP
    #   Keep only SNP-type variants. Indels, mixed sites, and
    #   symbolic alleles are discarded.
    #
    # --restrict-alleles-to BIALLELIC
    #   Keep only sites with exactly one alt allele.
    #   Multiallelic SNPs (e.g., A→C and A→T at the same site)
    #   are rare but problematic for downstream analyses.

    echo "SNP extraction complete: $(date)" | tee -a "$LOG"
fi

# Count variants for log
SNP_COUNT=$(gatk CountVariants -V "$SNP_VCF" 2>/dev/null | tail -1)
echo "Biallelic SNPs extracted: $SNP_COUNT" | tee -a "$LOG"

# ============================================================
# STEP 2 — VariantFiltration: tag sites that fail hard filters
# ============================================================
# Apply the GATK-recommended hard filter thresholds for SNPs.
# Sites that fail one or more filters are tagged — FILTER field
# is set to the filter name(s) they failed (e.g., "FS60;SOR3").
# Sites that pass all filters receive FILTER=PASS.
#
# IMPORTANT: VariantFiltration does NOT remove sites.
# The tagged VCF retains all variants so you can examine which
# filters fired most often before the final removal in Step 3.

echo "" | tee -a "$LOG"
echo "=== Step 2: Applying hard filters ===" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"

if [[ -f "$TAGGED_VCF" ]]; then
    echo "Tagged VCF already exists: $TAGGED_VCF — skipping Step 2." | tee -a "$LOG"
else
    gatk --java-options "-Xmx16g -Djava.io.tmpdir=$OUT/tmp" VariantFiltration \
        -R "$REF" \
        -V "$SNP_VCF" \
        -O "$TAGGED_VCF" \
        --filter-expression "QD < 2.0"            --filter-name "QD2" \
        --filter-expression "QUAL < 30.0"         --filter-name "QUAL30" \
        --filter-expression "FS > 60.0"           --filter-name "FS60" \
        --filter-expression "SOR > 3.0"           --filter-name "SOR3" \
        --filter-expression "MQ < 40.0"           --filter-name "MQ40" \
        --filter-expression "MQRankSum < -12.5"   --filter-name "MQRankSum-12.5" \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        2>&1 | tee -a "$LOG"
    #
    # Each --filter-expression / --filter-name pair defines one filter:
    #   --filter-expression : the condition that flags a site as FAILING
    #   --filter-name       : the label written to the FILTER field if it fails
    #
    # QD2         : QualByDepth < 2.0
    #               Variant quality is low relative to read depth.
    #               Likely a weak call with insufficient support.
    #
    # QUAL30      : QUAL < 30.0
    #               Phred-scaled confidence that the variant exists at this
    #               site. QUAL 30 = 1-in-1000 error probability. Low QUAL
    #               means GATK has weak overall confidence the variant is
    #               real, independent of read depth (unlike QD).
    #
    # FS60        : FisherStrand > 60.0
    #               Strong strand bias by Fisher's exact test.
    #               Variants supported by reads on only one strand
    #               are likely PCR or library prep artifacts.
    #
    # SOR3        : StrandOddsRatio > 3.0
    #               Strand bias by the more robust SOR metric.
    #               Applied alongside FS for complementary coverage.
    #
    # MQ40        : RMSMappingQuality < 40.0
    #               Reads at this site map poorly to the reference.
    #               Low mapping quality suggests a repetitive region
    #               where reads are not placed reliably.
    #
    # MQRankSum-12.5 : MappingQualityRankSum < -12.5
    #               Alt-allele reads have much lower mapping quality
    #               than ref-allele reads — the alt allele is supported
    #               only by poorly-mapping reads, which is an artifact
    #               signal rather than a real variant.
    #
    # ReadPosRankSum-8 : ReadPosRankSum < -8.0
    #               Alt allele tends to appear near the ends of reads,
    #               where sequencing errors are more common.

    echo "Hard filter tagging complete: $(date)" | tee -a "$LOG"
fi

# Summarise how many sites failed each filter
# Each tagged site has the filter name(s) it failed in the FILTER column (column 7).
# Sites failing multiple filters have semicolon-separated names (e.g., "FS60;SOR3").
# zcat decompresses the bgzipped VCF; grep -v "^#" skips header lines;
# awk counts rows where column 7 contains the filter name.
echo "" | tee -a "$LOG"
echo "Filter failure summary (sites failing each filter):" | tee -a "$LOG"
for filter in QD2 QUAL30 FS60 SOR3 MQ40 MQRankSum-12.5 ReadPosRankSum-8; do
    count=$(zcat "$TAGGED_VCF" | grep -v "^#" | awk -F'\t' -v f="$filter" '$7 ~ f {n++} END {print n+0}')
    echo "  $filter: $count sites" | tee -a "$LOG"
done

# ============================================================
# STEP 3 — SelectVariants: keep PASS sites only
# ============================================================
# Remove all sites that were tagged in Step 2, keeping only
# FILTER=PASS sites. This is the final filtered SNP dataset
# that will be passed to script 09 for population-level filtering.

echo "" | tee -a "$LOG"
echo "=== Step 3: Extracting PASS sites ===" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"

if [[ -f "$FILTERED_VCF" ]]; then
    echo "Filtered VCF already exists: $FILTERED_VCF — skipping Step 3." | tee -a "$LOG"
else
    gatk --java-options "-Xmx16g -Djava.io.tmpdir=$OUT/tmp" SelectVariants \
        -R "$REF" \
        -V "$TAGGED_VCF" \
        -O "$FILTERED_VCF" \
        --exclude-filtered \
        2>&1 | tee -a "$LOG"
    #
    # --exclude-filtered
    #   Remove any site where the FILTER field is not "PASS".
    #   This discards all sites tagged in Step 2, leaving only
    #   variants that passed all six hard filters.

    echo "PASS extraction complete: $(date)" | tee -a "$LOG"
fi

# Count final PASS variants
PASS_COUNT=$(gatk CountVariants -V "$FILTERED_VCF" 2>/dev/null | tail -1)
echo "PASS SNPs retained: $PASS_COUNT" | tee -a "$LOG"

# ============================================================
# SUMMARY
# ============================================================
echo "" | tee -a "$LOG"
echo "==================================================" | tee -a "$LOG"
echo "=== All steps complete ===" | tee -a "$LOG"
echo "Finished: $(date)" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "Biallelic SNPs before filtering : $SNP_COUNT" | tee -a "$LOG"
echo "PASS SNPs after hard filters    : $PASS_COUNT" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "SNPs only (pre-filter)  : $SNP_VCF" | tee -a "$LOG"
echo "Tagged VCF (for review) : $TAGGED_VCF" | tee -a "$LOG"
echo "Filtered VCF (PASS only): $FILTERED_VCF" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "Next step: script 09 — filter QC (09_filter_qc.Rmd); evaluate annotation distributions and Ti/Tv before population-level filtering (script 10)" | tee -a "$LOG"
