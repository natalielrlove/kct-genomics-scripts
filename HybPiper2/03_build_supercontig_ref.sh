#!/bin/bash
# ============================================================
# 03_build_supercontig_ref.sh — Build clean supercontig reference
# ============================================================
# Filters paralog-flagged loci and concatenates clean supercontig
# sequences into a single reference FASTA for downstream variant
# calling.
#
# WHAT IS A SUPERCONTIG REFERENCE?
# ---------------------------------
# hybpiper retrieve_sequences (run in 01_hybpiper_assemble.sh)
# produced one multi-sample FASTA per locus — e.g. 7271_supercontig.FNA
# contains the supercontig sequence (exon + flanking intron) for
# all 61 samples at locus 7271.
#
# For variant calling, we need a single representative sequence
# per locus as the mapping reference. Following Slimp et al. (2021,
# Applications in Plant Sciences), we select the LONGEST supercontig
# across all samples for each locus. The longest sequence has the
# most flanking intron recovered and gives all samples the best
# possible mapping target.
#
# This reference is much smaller than the whole Gymnocladus genome
# (~1.5–2 Mb across ~928 clean loci vs. 710 Mb for the full
# assembly). Variant calling against it focuses entirely on the
# target loci and runs much faster.
#
# PARALOG FILTERING
# -----------------
# paralog_report.tsv (from 02_paralog_retriever.sh) records the
# number of sequences recovered per sample × locus:
#   0 = locus not recovered in this sample
#   1 = one sequence (main target, no paralog)
#   2+ = multiple sequences — a paralog was assembled alongside
#        the main target, indicating a duplicated gene
#
# We use a conservative cutoff: exclude any locus where at least
# one sample has value > 1. This removes 77 / 1005 loci (7.7%),
# leaving 928 clean loci.
#
# OUTPUTS
# -------
#   paralog_loci.txt           — 77 excluded locus IDs
#   clean_loci.txt             — 928 retained locus IDs
#   gymno_supercontigs_clean.fasta — reference FASTA (one seq
#                                    per locus, 928 entries)
#
# Requirements: conda environment "hybpiper" must be active
#   source /home/nlove/miniconda3/etc/profile.d/conda.sh
#   conda activate hybpiper
# ============================================================

set -euo pipefail

# ============================================================
# --- Configuration -------------------------------------------
# ============================================================

# HybPiper output directory — contains supercontig FNA files
# and paralog_report.tsv from previous steps
OUT="/home/nlove/kct_genomics/output_files/hybpiper"

# Paralog report from 02_paralog_retriever.sh
PARALOG_REPORT="$OUT/paralog_report.tsv"

# Output reference FASTA
REF="$OUT/gymno_supercontigs_clean.fasta"

# ============================================================
# --- STEP 1: Identify paralog-flagged loci -------------------
# ============================================================
# Read paralog_report.tsv and flag any locus (column) where at
# least one sample (row) has a value > 1.
# Value > 1 means multiple sequences were recovered for that
# locus in that sample, indicating a paralog.

echo "============================================================"
echo "STEP 1: Identifying paralog-flagged loci"
echo "Paralog report: $PARALOG_REPORT"
echo "============================================================"

PARALOG_LOCI="$OUT/paralog_loci.txt"
CLEAN_LOCI="$OUT/clean_loci.txt"

awk '
NR == 1 {
    for (i = 2; i <= NF; i++) header[i] = $i
    ncols = NF
}
NR > 1 {
    for (i = 2; i <= NF; i++)
        if ($i > 1) paralog[i] = 1
}
END {
    for (i = 2; i <= ncols; i++) {
        if (i in paralog)
            print header[i] > "/dev/stderr"
        else
            print header[i]
    }
}
' "$PARALOG_REPORT" > "$CLEAN_LOCI" 2> "$PARALOG_LOCI"

N_PARALOG=$(wc -l < "$PARALOG_LOCI")
N_CLEAN=$(wc -l < "$CLEAN_LOCI")

echo "Paralog-flagged loci (excluded): $N_PARALOG"
echo "Clean loci (retained)          : $N_CLEAN"
echo "Written to:"
echo "  $PARALOG_LOCI"
echo "  $CLEAN_LOCI"

# ============================================================
# --- STEP 2: Build reference FASTA --------------------------
# ============================================================
# For each clean locus, select the LONGEST supercontig sequence
# across all samples and write it to the reference FASTA.
#
# This follows Slimp et al. (2021, Applications in Plant Sciences):
#   "reference sequences were constructed by selecting the longest
#    supercontig recovered for that species for each gene" (Fig. 2).
#
# Rationale: the longest supercontig has the most flanking intron
# sequence recovered, giving reads from all samples the best
# possible mapping target. Using the first sample's sequence
# would be arbitrary and potentially shorter.
#
# The supercontig FNA files produced by retrieve_sequences are
# named: {locus}_supercontig.FNA
# Each file is a multi-FASTA with one entry per sample.

echo ""
echo "============================================================"
echo "STEP 2: Building clean supercontig reference FASTA"
echo "(selecting longest supercontig per locus, per Slimp et al. 2021)"
echo "============================================================"

> "$REF"   # empty the output file before writing

N_WRITTEN=0
N_MISSING=0

while IFS= read -r locus; do
    fna="$OUT/${locus}_supercontig.FNA"

    if [[ ! -f "$fna" ]]; then
        echo "WARNING: missing FNA for locus $locus — skipping"
        (( N_MISSING++ )) || true
        continue
    fi

    # Extract the longest sequence from this multi-FASTA.
    # awk accumulates each sequence, tracks the max length,
    # and prints only the header + sequence of the longest entry.
    awk '
    /^>/ {
        if (seq != "") {
            seqlen = length(seq)
            if (seqlen > maxlen) {
                maxlen = seqlen
                maxheader = header
                maxseq = seq
            }
        }
        header = $0
        seq = ""
        next
    }
    { seq = seq $0 }
    END {
        if (length(seq) > maxlen) {
            maxheader = header
            maxseq = seq
        }
        print maxheader
        print maxseq
    }
    ' "$fna" >> "$REF"

    (( N_WRITTEN++ )) || true

done < "$CLEAN_LOCI"

echo "Loci written to reference : $N_WRITTEN"
if [[ "$N_MISSING" -gt 0 ]]; then
    echo "WARNING: $N_MISSING loci had no FNA file and were skipped"
fi

# ============================================================
# --- STEP 3: Verify and index the reference -----------------
# ============================================================

echo ""
echo "============================================================"
echo "STEP 3: Indexing reference with BWA and samtools"
echo "============================================================"

# Count sequences in the reference
N_SEQS=$(grep -c "^>" "$REF")
echo "Sequences in reference: $N_SEQS"

# BWA index — required for read mapping in 04_variantcall.sh
bwa index "$REF"

# samtools fai index — required by GATK HaplotypeCaller
samtools faidx "$REF"

# GATK sequence dictionary — required by GATK HaplotypeCaller.
# Creates a .dict file listing contig names and lengths.
# GATK requires all three: .fai, .dict, and BWA index.
gatk CreateSequenceDictionary -R "$REF"

echo ""
echo "============================================================"
echo "Reference build complete."
echo "Key outputs:"
echo "  Clean locus list  : $CLEAN_LOCI"
echo "  Paralog locus list : $PARALOG_LOCI"
echo "  Reference FASTA   : $REF"
echo "  BWA index         : ${REF}.{amb,ann,bwt,pac,sa}"
echo "  samtools index    : ${REF}.fai"
echo "  GATK dictionary   : ${REF%.fasta}.dict"
echo "============================================================"
