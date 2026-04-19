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
| `gatk` | `environment_gatk.yml` | scripts 05–08 (gatk4, samtools) |

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
| Chromosomes (`chr*`) | 14 | 659,219,136 | ~93% |
| Unplaced scaffolds (`scaffold_*`) | 972 | 51,323,761 | ~7% |
| **Total** | **986** | **710,542,897** | |

Scaffold sizes range from ~660 kb (largest) to ~25 kb (smallest), averaging ~57 kb.
For probe design: consider excluding unplaced scaffolds entirely or setting a minimum
size cutoff (e.g. >100 kb) to avoid repetitive/low-complexity regions.

## Pipeline run log

| Script | Start | End | Wall time | Notes |
|--------|-------|-----|-----------|-------|
| `06_haplotypecaller.sh` | 2026-03-18 ~16:00 PDT | 2026-03-23 ~16:25 PDT | ~120 h | 61 samples, 10 parallel jobs × 3 threads, Hap 1 reference (986 scaffolds); all 61 GVCFs + .tbi produced, no errors |
| `07_genomicsdb_import.sh` | 2026-03-26 | 2026-03-28 ~16:33 UTC | ~41 h (2,471.61 min) | GenomicsDBImport + GenotypeGVCFs across 986 scaffolds; outputs: genomicsdb workspace + `gymno_hap1.raw.vcf.gz` |
| `08_filter_snps.sh` | 2026-04-01 03:19 UTC | 2026-04-01 04:32 UTC | ~1h 13min | 16,021,590 biallelic SNPs → 13,672,054 PASS (14.7% removed); outputs in `gatk_filter_hap1/` |

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
