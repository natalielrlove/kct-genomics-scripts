# kct-genomics-scripts

Analysis scripts for *Gymnocladus dioicus* population genomics.

## Repository structure

```
scripts/
├── 00_file_organization/    # File organization and setup
├── 01_probe_validation/     # In silico probe validation (hyb-seq)
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
