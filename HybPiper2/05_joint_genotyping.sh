#!/bin/bash
# ============================================================
# 05_joint_genotyping.sh — Joint genotype + hard-filter HybPiper gVCFs
# ============================================================
# Takes the 61 per-sample gVCFs produced by 04_variantcall.sh
# and produces a filtered multi-sample SNP VCF through five steps:
#
#   Step 1 — Generate sample map
#   Step 2 — Generate interval list from supercontig reference index
#   Step 3 — GenomicsDBImport: load all 61 gVCFs into a database
#   Step 4 — GenotypeGVCFs: joint-genotype across all 61 samples
#   Step 5 — SelectVariants: extract biallelic SNPs only
#   Step 6 — VariantFiltration: tag SNPs that fail hard filters
#   Step 7 — SelectVariants: keep PASS sites only
#
# WHY COMBINE GENOTYPING AND FILTERING HERE?
# ------------------------------------------
# In the WGS pipeline these were two scripts (07 + 08) because
# the genome-wide GenomicsDBImport + GenotypeGVCFs took ~41 h,
# after which we wanted a natural stopping point to inspect the
# raw VCF before filtering. The supercontig reference is ~375×
# smaller (~1.9 Mb vs 710 Mb), so the full pipeline runs in
# hours rather than days — combining the steps is practical and
# matches the Slimp et al. (2021) variantcall.sh structure.
#
# RUNTIME ESTIMATE (supercontig reference ~1.9 Mb / 928 loci)
# ------------------------------------------------------------
# GenomicsDBImport : ~30–60 min (61 gVCFs, small reference)
# GenotypeGVCFs   : ~15–30 min
# Filtering steps : ~5 min each
# Total           : ~1.5–2 h
#
# OUTPUT
# ------
# GenomicsDB workspace:
#   output_files/hybpiper_joint/genomicsdb/
#
# Raw joint VCF:
#   output_files/hybpiper_joint/gymno_sc.raw.vcf.gz
#
# Biallelic SNPs (pre-filter, for inspection):
#   output_files/hybpiper_joint/gymno_sc.snps.vcf.gz
#
# SNPs with filter tags (for review):
#   output_files/hybpiper_joint/gymno_sc.snps.tagged.vcf.gz
#
# Final PASS SNPs (input for 06_plink.sh):
#   output_files/hybpiper_joint/gymno_sc.snps.filtered.vcf.gz
#
# RESUME BEHAVIOR
# ---------------
# GenomicsDBImport exits with a clear message if the workspace
# already exists (cannot be overwritten). All other steps are
# skipped if their output file already exists.
# To re-run from scratch:
#   rm -rf output_files/hybpiper_joint/genomicsdb/
#
# Requirements: conda environment "gatk" must be active
#   source /home/nlove/miniconda3/etc/profile.d/conda.sh
#   conda activate gatk
# ============================================================

set -euo pipefail

# ============================================================
# --- Paths ---------------------------------------------------
# ============================================================

# Supercontig reference built by 03_build_supercontig_ref.sh
REF="/home/nlove/kct_genomics/output_files/hybpiper/gymno_supercontigs_clean.fasta"

# Input: per-sample gVCFs from 04_variantcall.sh
# These landed on /data after the pipeline completed.
GVCF_DIR="/data/labs/Love/Shared/output_files/hybpiper_variantcall/gvcf"

# Output: joint genotyping and filtering results
# Always write to /home (local disk, faster I/O). Move to /data only after pipeline finishes.
OUT="/home/nlove/kct_genomics/output_files/hybpiper_joint"

# --- Output file names ---------------------------------------
WORKSPACE="$OUT/genomicsdb"
SAMPLE_MAP="$OUT/tmp/sample_map.tsv"
INTERVALS="$OUT/tmp/all_loci.list"
RAW_VCF="$OUT/gymno_sc.raw.vcf.gz"
SNP_VCF="$OUT/gymno_sc.snps.vcf.gz"
TAGGED_VCF="$OUT/gymno_sc.snps.tagged.vcf.gz"
FILTERED_VCF="$OUT/gymno_sc.snps.filtered.vcf.gz"
LOG="$OUT/logs/05_joint_genotyping.log"

mkdir -p "$OUT"/{logs,tmp}

echo "============================================================" | tee "$LOG"
echo "05_joint_genotyping.sh" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"
echo "Reference  : $REF" | tee -a "$LOG"
echo "gVCF input : $GVCF_DIR" | tee -a "$LOG"
echo "Output dir : $OUT" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"

# ============================================================
# STEP 1 — Generate sample map
# ============================================================
# GenomicsDBImport requires a tab-separated file mapping each
# sample name to its gVCF path:
#   Gymno_2013_109 <TAB> /path/to/Gymno_2013_109.sc.g.vcf.gz
#
# Sample name is the .sc.g.vcf.gz filename stripped of the
# directory and extension, matching the SM tag used in script 04.

echo "" | tee -a "$LOG"
echo "=== Step 1: Generating sample map ===" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"

for gvcf in "$GVCF_DIR"/*.sc.g.vcf.gz; do
    sample=$(basename "$gvcf" .sc.g.vcf.gz)
    printf "%s\t%s\n" "$sample" "$gvcf"
done > "$SAMPLE_MAP"

echo "Sample map written to : $SAMPLE_MAP" | tee -a "$LOG"
echo "Samples in map        : $(wc -l < "$SAMPLE_MAP")" | tee -a "$LOG"

# ============================================================
# STEP 2 — Generate interval list
# ============================================================
# GenomicsDBImport requires -L (intervals). We include all loci
# by extracting scaffold names from the reference .fai index.
# The supercontig reference has 928 loci (one contig per locus).

echo "" | tee -a "$LOG"
echo "=== Step 2: Generating interval list ===" | tee -a "$LOG"

awk '{print $1}' "${REF}.fai" > "$INTERVALS"

echo "Interval list written to : $INTERVALS" | tee -a "$LOG"
echo "Loci in interval list    : $(wc -l < "$INTERVALS")" | tee -a "$LOG"

# ============================================================
# STEP 3 — GenomicsDBImport
# ============================================================
# Loads all 61 per-sample gVCFs into a GenomicsDB data store.
# --merge-input-intervals combines all 928 loci into one database
# partition rather than creating a separate partition per locus
# (recommended for fragmented assemblies like this supercontig set).
# --reader-threads 4 reads input gVCFs in parallel before the
# single-threaded merge.

echo "" | tee -a "$LOG"
echo "=== Step 3: GenomicsDBImport ===" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"

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

echo "" | tee -a "$LOG"
echo "GenomicsDBImport complete: $(date)" | tee -a "$LOG"
echo "Database written to: $WORKSPACE" | tee -a "$LOG"

# ============================================================
# STEP 4 — GenotypeGVCFs
# ============================================================
# Queries the GenomicsDB database and performs joint genotyping
# across all 61 samples. Outputs a raw multi-sample VCF with
# one row per variant site and one column per sample.
# The "gendb://" prefix tells GATK to read from a GenomicsDB
# store rather than a plain VCF path.

echo "" | tee -a "$LOG"
echo "=== Step 4: GenotypeGVCFs ===" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"

if [[ -f "$RAW_VCF" ]]; then
    echo "Raw VCF already exists: $RAW_VCF — skipping Step 4." | tee -a "$LOG"
else
    gatk --java-options "-Xmx64g -Djava.io.tmpdir=$OUT/tmp" GenotypeGVCFs \
        -R "$REF" \
        -V "gendb://$WORKSPACE" \
        -O "$RAW_VCF" \
        --tmp-dir "$OUT/tmp" \
        2>&1 | tee -a "$LOG"

    echo "GenotypeGVCFs complete: $(date)" | tee -a "$LOG"
fi

echo "Raw joint VCF: $RAW_VCF" | tee -a "$LOG"

# ============================================================
# STEP 5 — SelectVariants: extract biallelic SNPs
# ============================================================
# Extract SNP-type variants only and restrict to biallelic sites.
# Multiallelic SNPs are rare but problematic for PLINK and most
# population genomics tools, which assume exactly two alleles.

echo "" | tee -a "$LOG"
echo "=== Step 5: Extracting biallelic SNPs ===" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"

if [[ -f "$SNP_VCF" ]]; then
    echo "SNP VCF already exists: $SNP_VCF — skipping Step 5." | tee -a "$LOG"
else
    gatk --java-options "-Xmx16g -Djava.io.tmpdir=$OUT/tmp" SelectVariants \
        -R "$REF" \
        -V "$RAW_VCF" \
        -O "$SNP_VCF" \
        --select-type-to-include SNP \
        --restrict-alleles-to BIALLELIC \
        2>&1 | tee -a "$LOG"

    echo "SNP extraction complete: $(date)" | tee -a "$LOG"
fi

SNP_COUNT=$(gatk CountVariants -V "$SNP_VCF" 2>/dev/null | tail -1)
echo "Biallelic SNPs extracted: $SNP_COUNT" | tee -a "$LOG"

# ============================================================
# STEP 6 — VariantFiltration: tag sites failing hard filters
# ============================================================
# Apply GATK-recommended hard filter thresholds for SNPs.
# Sites that fail are tagged in the FILTER field (e.g., "FS60");
# sites that pass all filters receive FILTER=PASS.
# This step does NOT remove sites — tagging first lets you inspect
# which filters fire before committing to removal in Step 7.
#
# THRESHOLDS (same as WGS script 08, GATK recommended for SNPs)
# These are starting points — check annotation distributions in
# the script 09 QC notebook and adjust if a filter removes >50%
# of sites or if annotation peaks look unusual for this dataset.
#
#   QD   < 2.0    : Variant quality low relative to read depth
#   QUAL < 30.0   : Overall call confidence < 1-in-1000
#   FS   > 60.0   : Strong strand bias (Fisher's exact test)
#   SOR  > 3.0    : Strand bias (odds ratio; complements FS)
#   MQ   < 40.0   : Poor read mapping quality at this site
#   MQRankSum < -12.5 : Alt allele reads map worse than ref reads
#   ReadPosRankSum < -8.0 : Alt allele enriched at read ends

echo "" | tee -a "$LOG"
echo "=== Step 6: Applying hard filters ===" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"

if [[ -f "$TAGGED_VCF" ]]; then
    echo "Tagged VCF already exists: $TAGGED_VCF — skipping Step 6." | tee -a "$LOG"
else
    gatk --java-options "-Xmx16g -Djava.io.tmpdir=$OUT/tmp" VariantFiltration \
        -R "$REF" \
        -V "$SNP_VCF" \
        -O "$TAGGED_VCF" \
        --filter-expression "QD < 2.0"              --filter-name "QD2" \
        --filter-expression "QUAL < 30.0"           --filter-name "QUAL30" \
        --filter-expression "FS > 60.0"             --filter-name "FS60" \
        --filter-expression "SOR > 3.0"             --filter-name "SOR3" \
        --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
        --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        2>&1 | tee -a "$LOG"

    echo "Hard filter tagging complete: $(date)" | tee -a "$LOG"
fi

# Report how many sites failed each filter
echo "" | tee -a "$LOG"
echo "Filter failure summary (sites failing each filter):" | tee -a "$LOG"
for filter in QD2 QUAL30 FS60 SOR3 MQ40 MQRankSum-12.5 ReadPosRankSum-8; do
    count=$(zcat "$TAGGED_VCF" | grep -v "^#" | awk -F'\t' -v f="$filter" '$7 ~ f {n++} END {print n+0}')
    echo "  $filter: $count sites" | tee -a "$LOG"
done

# ============================================================
# STEP 7 — SelectVariants: keep PASS sites only
# ============================================================
# Remove all tagged sites, retaining only FILTER=PASS variants.
# This is the final filtered SNP set passed to 06_plink.sh.

echo "" | tee -a "$LOG"
echo "=== Step 7: Extracting PASS sites ===" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"

if [[ -f "$FILTERED_VCF" ]]; then
    echo "Filtered VCF already exists: $FILTERED_VCF — skipping Step 7." | tee -a "$LOG"
else
    gatk --java-options "-Xmx16g -Djava.io.tmpdir=$OUT/tmp" SelectVariants \
        -R "$REF" \
        -V "$TAGGED_VCF" \
        -O "$FILTERED_VCF" \
        --exclude-filtered \
        2>&1 | tee -a "$LOG"

    echo "PASS extraction complete: $(date)" | tee -a "$LOG"
fi

PASS_COUNT=$(gatk CountVariants -V "$FILTERED_VCF" 2>/dev/null | tail -1)
echo "PASS SNPs retained: $PASS_COUNT" | tee -a "$LOG"

# ============================================================
# SUMMARY
# ============================================================
echo "" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
echo "=== All steps complete ===" | tee -a "$LOG"
echo "Finished: $(date)" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "Biallelic SNPs before filtering : $SNP_COUNT" | tee -a "$LOG"
echo "PASS SNPs after hard filters    : $PASS_COUNT" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "GenomicsDB workspace    : $WORKSPACE" | tee -a "$LOG"
echo "Raw joint VCF           : $RAW_VCF" | tee -a "$LOG"
echo "SNPs only (pre-filter)  : $SNP_VCF" | tee -a "$LOG"
echo "Tagged VCF (for review) : $TAGGED_VCF" | tee -a "$LOG"
echo "Filtered VCF (PASS only): $FILTERED_VCF" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "Next step: 06_plink.sh — LD pruning, PCA, diversity stats" | tee -a "$LOG"
echo "============================================================" | tee -a "$LOG"
