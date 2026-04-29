# kct-genomics-scripts

Analysis scripts for *Gymnocladus dioicus* population genomics.

## Repository structure

```
kct-genomics-scripts/
├── 00_file_organization/     # File organization and setup
├── 01_QA_QC/                 # Read quality assessment and trimming  [env: qc]
│   ├── 01_fastqc.sh          # Per-sample quality reports (FastQC)
│   ├── 02_multiqc.sh         # Aggregate FastQC reports (MultiQC)
│   └── 03_fastp.sh           # Adapter trimming and filtering (fastp)
├── 02_mapping/               # Alignment to reference genome          [env: mapping]
│   └── 04_bwamem.sh          # Paired-end alignment with bwa-mem2 (2 parallel jobs, 15 threads each)
├── 03_calling_snps/          # GATK variant calling                   [env: gatk]
│   ├── 05_gatk_markdups.sh   # Add read groups + mark PCR duplicates (5 parallel jobs)
│   ├── 06_haplotypecaller.sh # Per-sample variant calling in gVCF mode (10 parallel jobs, 3 threads each)
│   ├── 07_genomicsdb_import.sh # Joint genotyping: GenomicsDBImport + GenotypeGVCFs
│   └── 08_filter_snps.sh     # GATK hard filters: extract biallelic SNPs → tag → keep PASS only
├── 04_filter_qc/             # Hard filter QC and threshold evaluation  [env: R]
│   └── 09_filter_qc.Rmd      # Annotation distributions, Ti/Tv ratio, threshold review
├── HybPiper2/                # Target locus recovery from WGS data   [env: hybpiper]
│   ├── 00_prep_targets.sh    # One-time reformat of Crameri target FASTA for HybPiper2 (run locally)
│   └── 01_hybpiper_assemble.sh # Full HybPiper2 pipeline: assemble → stats → heatmap → retrieve
├── Fabaceae_probes_targets/  # Target sequences for HybPiper2 (committed to repo)
│   ├── Fabaceae_iter2_1005reg.fasta          # Original Crameri et al. sequences (1,005 loci)
│   └── Fabaceae_iter2_1005reg_hybpiper.fasta # HybPiper-ready version (reformatted headers)
├── environment_qc.yml        # Conda environment for scripts 01-03 (fastqc, fastp, multiqc)
├── environment_mapping.yml   # Conda environment for script 04 (bwa-mem2, samtools, parallel)
├── environment_gatk.yml      # Conda environment for scripts 05-07 (gatk4, samtools)
└── README.md
```

## Data structure (on server, not tracked by git)

```
~/kct_genomics/
├── scripts/           ← this repo
├── input_files/       ← raw data (genomes, FASTQs, probe FASTA)
└── output_files/      ← pipeline outputs, organized by step
```

## Server
Scripts are run on the Curie server. Clone this repo to `~/kct_genomics/scripts/`.
Push changes from local machine, pull on server to run.

## Conda environments

Four conda environments are used across the pipeline. The yml files in this repo
allow anyone to recreate them exactly.

| Environment | yml file | Used by |
|---|---|---|
| `qc` | `environment_qc.yml` | scripts 01–03 (fastqc, fastp, multiqc) |
| `mapping` | `environment_mapping.yml` | script 04 (bwa-mem2, samtools, parallel) |
| `gatk` | `environment_gatk.yml` | scripts 05–08 (gatk4, samtools) |
| `hybpiper` | `environment_hybpiper.yml` | HybPiper2 scripts (hybpiper, bwa, spades, parallel) |

**To recreate an environment from scratch:**
```bash
conda env create -f environment_gatk.yml
```

**To activate before running scripts:**
```bash
source /home/nlove/miniconda3/etc/profile.d/conda.sh
conda activate gatk   # or qc, or mapping
```

Each script's header specifies which environment it requires.

## HybPiper2 — Target locus recovery

HybPiper2 is used to assess how well the Crameri et al. (2022) **Fabaceae1005** probe set recovers target loci from *Gymnocladus dioicus* whole-genome sequencing data. The approach treats 30× WGS reads as if they came from a hybridization capture experiment (virtual capture), providing a recovery estimate before committing to wet-lab probe synthesis.

**Probe set:** 1,005 conserved nuclear loci derived from the *Cajanus cajan* (pigeon pea) reference genome, validated across all six Fabaceae subfamilies including Caesalpinioideae (*Gymnocladus*).

### Scripts

| Script | Purpose |
|--------|---------|
| `HybPiper2/00_prep_targets.sh` | One-time local reformat of Crameri target FASTA for HybPiper2 compatibility (run on Mac before first use) |
| `HybPiper2/01_hybpiper_assemble.sh` | Full pipeline: assemble all samples → stats → recovery heatmap → retrieve sequences |

### Target file format — CRITICAL

HybPiper2 requires target FASTA headers in the format `>taxon-gene`, where the **part after the last dash** becomes the locus directory name. With 1,005 loci this must be:

```
>Cajanus_cajan-7271    ← correct: gene=7271 → 1005 unique locus directories
>Cajanus_cajan-243
```

**Not** the reverse (`>7271-Cajanus_cajan`), which causes HybPiper to use `Cajanus_cajan` as the gene name for all sequences → one directory instead of 1005.

`00_prep_targets.sh` produces the correct format. The reformatted file (`Fabaceae_iter2_1005reg_hybpiper.fasta`) is committed to this repo and ready to use.

### Expected warnings during assembly

These warnings are normal and do not indicate errors:

- **Stop codons in translated sequence** — expected for genomic/non-CDS targets that include introns; these are conserved windows, not annotated exons
- **Paralog warnings** — normal for a diploid species; HybPiper selects the most likely ortholog
- **Non-consecutive Exonerate hits** — rare; affects supercontig stitching only, not the primary exon sequence

Confirmation of a successful run: `"Generated sequences from 1005 genes!"` in the assembly log.

### Key outputs

| File/Directory | Description |
|---|---|
| `recovery_heatmap.png` | Heatmap of % target recovered per sample × locus — primary diagnostic |
| `seq_lengths.tsv` | Numeric recovery table (input to heatmap) |
| `*.FNA` | One multi-FASTA per locus across all samples (analysis-ready) |
| `*_supercontig.fasta` | Exon + flanking intron sequences (more variable; better for within-species analyses) |

## Reference genome assembly stats

Reference: `Gymnocladus_dioicus_M_hap1.fa` (Haplotype 1, 986 sequences)

| Sequence type | Count | Total bp | % of genome |
|---|---|---|---|
| Chromosomes (`chr*`) | 14 | 659,219,136 | ~93% |
| Unplaced scaffolds (`scaffold_*`) | 972 | 51,323,761 | ~7% |
| **Total** | **986** | **710,542,897** | |

Scaffold sizes range from ~660 kb (largest) to ~25 kb (smallest), averaging ~57 kb.
For probe design: consider excluding unplaced scaffolds entirely or setting a minimum
size cutoff (e.g. >100 kb) to avoid repetitive/low-complexity regions.

## Pipeline run log

| Script | Start | End | Wall time | Notes |
|--------|-------|-----|-----------|-------|
| `01_hybpiper_assemble.sh` | 2026-04-23 ~02:31 UTC | in progress | — | 61 samples, 8 parallel jobs × 2 threads; correct target format `>Cajanus_cajan-7271`; numbered locus directories confirmed per sample |
| `05_gatk_markdups.sh` *(mq30 re-run)* | 2026-04-19 | 2026-04-20 | ~191 min (~3h 11min) | 10 parallel jobs; input: MQ30-filtered BAMs from `bam_mq30/` |
| `06_haplotypecaller.sh` *(mq30 re-run)* | ~2026-04-20 UTC | ~2026-04-22 22:38 UTC | ~4,299 min (~71.7h) | 61 samples, 10 parallel jobs × 3 threads, Hap 1 reference (986 scaffolds); faster than original run (~120h) due to MQ30 pre-filter reducing input reads |
| `07_genomicsdb_import.sh` *(mq30 re-run)* | 2026-04-22 22:38 UTC | 2026-04-25 06:44 UTC | ~3,365 min (~56.1h) | GenomicsDBImport: 1,656 min; GenotypeGVCFs: 1,708 min; outputs: genomicsdb workspace + `gymno_hap1.raw.vcf.gz`; original run was ~41h (2,471 min) |
| `08_filter_snps.sh` *(mq30 re-run)* | 2026-04-25 06:44 UTC | 2026-04-25 07:26 UTC | ~42 min | 15,330,896 biallelic SNPs → 13,710,712 PASS (10.6% removed); original run: 16,021,590 → 13,672,054 PASS (14.7% removed) |

**Script 08 filter failure summary** (sites failing each filter; sites can fail multiple filters):

| Filter | Sites failed |
|--------|-------------|
| QD2 (QualByDepth < 2.0) | 1,363,777 |
| QUAL30 (QUAL < 30.0) | 0 |
| FS60 (FisherStrand > 60.0) | 75,218 |
| SOR3 (StrandOddsRatio > 3.0) | 518,874 |
| MQ40 (RMSMappingQuality < 40.0) | 801,726 |
| MQRankSum-12.5 | 305 |
| ReadPosRankSum-8 | 4 |

**Script 08b missingness check** (run 2026-04-19, original pipeline, PASS SNPs only, chr01–chr14):

| Threshold | chr01 sites kept | chr01 total | chr01 % | chr01–14 sites kept | chr01–14 total | chr01–14 % |
|-----------|-----------------|-------------|---------|---------------------|----------------|------------|
| F_MISSING ≤ 0.5 | 1,195,784 | 1,212,145 | 98.6% | 13,462,101 | 13,662,332 | 98.5% |
| F_MISSING ≤ 0.3 | 1,176,855 | 1,212,145 | 97.0% | 13,168,342 | 13,662,332 | 96.3% |
| F_MISSING ≤ 0.1 | 1,108,361 | 1,212,145 | 91.4% | 12,077,197 | 13,662,332 | 88.3% |

Very low missingness overall — consistent with 30× WGS coverage. Planned threshold for script 10: F_MISSING ≤ 0.1.
