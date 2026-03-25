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
│   └── 07_genomicsdb_import.sh # Joint genotyping: GenomicsDBImport + GenotypeGVCFs
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

Three conda environments are used across the pipeline. The yml files in this repo
allow anyone to recreate them exactly.

| Environment | yml file | Used by |
|---|---|---|
| `qc` | `environment_qc.yml` | scripts 01–03 (fastqc, fastp, multiqc) |
| `mapping` | `environment_mapping.yml` | script 04 (bwa-mem2, samtools, parallel) |
| `gatk` | `environment_gatk.yml` | scripts 05–07 (gatk4, samtools) |

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

## Reference genome assembly stats

Reference: `Gymnocladus_dioicus_M_hap1.fa` (Haplotype 1, 986 sequences)

| Sequence type | Count | Total bp | % of genome |
|---|---|---|---|
| Chromosomes (`chr*`) | 86 | 659,219,136 | ~93% |
| Unplaced scaffolds (`scaffold_*`) | 900 | 51,323,761 | ~7% |
| **Total** | **986** | **710,542,897** | |

Scaffold sizes range from ~660 kb (largest) to ~25 kb (smallest), averaging ~57 kb.
For probe design: consider excluding unplaced scaffolds entirely or setting a minimum
size cutoff (e.g. >100 kb) to avoid repetitive/low-complexity regions.

## Pipeline run log

| Script | Start | End | Wall time | Notes |
|--------|-------|-----|-----------|-------|
| `06_haplotypecaller.sh` | 2026-03-18 ~16:00 PDT | 2026-03-23 ~16:25 PDT | ~120 h | 61 samples, 10 parallel jobs × 3 threads, Hap 1 reference (986 scaffolds); all 61 GVCFs + .tbi produced, no errors |
