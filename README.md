# kct-genomics-scripts

Analysis scripts for *Gymnocladus dioicus* population genomics.

## Repository structure

```
kct-genomics-scripts/
├── 00_file_organization/    # File organization and setup
├── 01_QA_QC/                # Read quality assessment and trimming
│   ├── 01_fastqc.sh         # Per-sample quality reports (FastQC)
│   ├── 02_multiqc.sh        # Aggregate FastQC reports (MultiQC)
│   └── 03_fastp.sh          # Adapter trimming and filtering (fastp)
├── 02_mapping/              # Alignment to reference genome
│   └── 04_bwamem.sh         # Paired-end alignment with bwa-mem2 (2 parallel jobs, 15 threads each)
├── 03_calling_snps/         # GATK variant calling
│   ├── 05_gatk_markdups.sh  # Add read groups + mark PCR duplicates (5 parallel jobs)
│   └── 06_haplotypecaller.sh # Per-sample variant calling in gVCF mode (10 parallel jobs, 3 threads each)
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

## Pipeline run log

| Script | Start | End | Wall time | Notes |
|--------|-------|-----|-----------|-------|
| `06_haplotypecaller.sh` | 2026-03-18 ~16:00 PDT | 2026-03-23 ~16:25 PDT | ~120 h | 61 samples, 10 parallel jobs × 3 threads, Hap 1 reference (986 scaffolds); all 61 GVCFs + .tbi produced, no errors |
