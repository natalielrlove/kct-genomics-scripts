#!/bin/bash
# ============================================================
# 00_prep_targets.sh — Reformat target FASTA for HybPiper2
# ============================================================
# Run this ONCE locally (on your Mac) before pushing to the
# server. It produces the HybPiper-ready target file from the
# Crameri et al. reference sequences.
#
# WHAT ARE THESE SEQUENCES?
# -------------------------
# The target file (Fabaceae_iter2_1005reg.fasta) contains 1,005
# reference sequences extracted from the Cajanus cajan (pigeon
# pea) reference genome. Cajanus cajan is a domesticated legume
# with one of the best-annotated genomes in Fabaceae, so Crameri
# used it as the source genome for designing the probes.
#
# Each sequence represents a conserved genomic region that was:
#   1. Present and single-copy across multiple legume genomes
#   2. Successfully captured in test hybridization experiments
#   3. Retained after two rounds of CaptureAL filtering
#      (iter1 = 2,468 regions after first pass;
#       iter2 = 1,005 regions after second, stricter pass)
#
# These are NOT necessarily annotated gene CDS sequences —
# they are conserved genomic windows (mean ~1,210 bp) that
# were found to work well for target capture across the family.
#
# WHAT IS NC_033804.1?
# --------------------
# This is a GenBank/RefSeq accession number — essentially a
# library call number for a specific chromosome sequence in the
# database. NC_033804.1 = chromosome 1 (LG01) of the Cajanus
# cajan reference genome. The coordinates (e.g., 1244349–1244570)
# are the position of the target region on that chromosome.
# Unanchored scaffolds that couldn't be placed on a chromosome
# are listed as NA_LG_Scaffold... instead of a NC_ accession.
#
# WHAT IS THE "GENE ID"?
# ----------------------
# HybPiper needs a short, stable identifier for each locus.
# The Crameri probe file uses IDs like "dalb_7271" — "dalb"
# was a development prefix (referencing Dalbergia, the focal
# genus), and 7271 is the target region number.
# The target file encodes this as "_ID_7271" at the end of
# each header. We extract just the number and use it as the
# HybPiper gene ID.
#
# WHY "Cajanus_cajan" AS THE TAXON NAME?
# ---------------------------------------
# HybPiper requires headers in the format >GENEID-taxon.
# The taxon label tells HybPiper which organism the reference
# sequence came from. Here it's Cajanus_cajan because these
# sequences were pulled from the pigeon pea genome.
# HybPiper uses this to handle paralogs — if you provide the
# same gene from multiple taxa, it can pick the best match.
# With only one taxon here, it functions as a simple reference.
# Your assembled Gymnocladus sequences will map against these
# Cajanus references and return Gymnocladus sequence — the
# taxon label in the target file does not affect your output.
#
# ORIGINAL HEADER FORMAT (from Crameri Dryad repository):
#   >NC_033804.1_LG_01_1244349_1244570_ID_7271
#   >NC_033804.1_LG_01_13802336_13802605_ID_7901.merged
#   >NA_LG_Scaffold098279_1_187_ID_243
#
# HYBPIPER-READY FORMAT (produced by this script):
#   >7271-Cajanus_cajan
#   >7901-Cajanus_cajan
#   >243-Cajanus_cajan
#
# IMPORTANT — FORMAT RULE:
#   HybPiper uses the part BEFORE the last dash as the locus/gene
#   identifier (what becomes the output directory name and locus ID
#   in all downstream stats). The part AFTER the last dash is the
#   taxon label. So gene ID MUST come first: >GENEID-TAXON.
#   Putting the taxon first (>Cajanus_cajan-7271) causes all 1005
#   sequences to be treated as ONE gene called "Cajanus_cajan".
#
# The ".merged" suffix (present on regions where two adjacent
# target windows were merged into one) is stripped — it is
# metadata about how the reference was assembled, not part of
# the locus identity.
# ============================================================

set -euo pipefail

# Input: Crameri iter2 reference sequences (as downloaded from Dryad)
IN="Fabaceae_probes_targets/Fabaceae_iter2_1005reg.fasta"

# Output: HybPiper-ready version with reformatted headers
OUT="Fabaceae_probes_targets/Fabaceae_iter2_1005reg_hybpiper.fasta"

echo "Reformatting headers for HybPiper2..."
echo "Input : $IN"
echo "Output: $OUT"

# --- Reformat headers ----------------------------------------
# For each header line, extract the ID number and append taxon.
#
# The sed pattern works as follows:
#   >.*_ID_   match the leading part of every header up to "_ID_"
#   \([0-9]*\) capture the numeric ID (e.g. 7271)
#   .*         match and discard everything after (including .merged)
#   >\1-Cajanus_cajan   replace with clean HybPiper format
#
# Sequence lines (not starting with >) are passed through unchanged.

sed 's/>.*_ID_\([0-9]*\).*/>\1-Cajanus_cajan/' "$IN" > "$OUT"

# --- Verify the output ---------------------------------------
N_IN=$(grep -c "^>" "$IN")
N_OUT=$(grep -c "^>" "$OUT")

echo ""
echo "Input sequences  : $N_IN"
echo "Output sequences : $N_OUT"

if [[ "$N_IN" -ne "$N_OUT" ]]; then
    echo "ERROR: sequence count mismatch — check the sed pattern"
    exit 1
fi

echo ""
echo "First 5 reformatted headers:"
grep "^>" "$OUT" | head -5

echo ""
echo "Done. Use $OUT as the -t_dna target file in 01_hybpiper_assemble.sh"
