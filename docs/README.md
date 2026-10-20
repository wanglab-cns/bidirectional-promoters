# Bidirectional Promoters in Pediatric Brain Tumors

**Authors:** [Amirfaham Fallahpour](https://github.com/Amirfaham1), [Farnoosh Abbas Aghababazadeh](https://github.com/RibaA), Kevin Wang

**Contact:** [bhklab.amirfahamfallahpour@gmail.com](mailto:bhklab.amirfahamfallahpour@gmail.com)

**Description:** A pan-brain screen of bidirectional promoters in pediatric brain tumors: which shared promoters are altered (methylation, copy number, mutation) and how each change affects both partner genes, including microRNAs.

--------------------------------------

[![pixi-badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json&style=flat-square)](https://github.com/prefix-dev/pixi)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json&style=flat-square)](https://github.com/astral-sh/ruff)
[![Built with Material for MkDocs](https://img.shields.io/badge/mkdocs--material-gray?logo=materialformkdocs&style=flat-square)](https://github.com/squidfunk/mkdocs-material)

![GitHub last commit](https://img.shields.io/github/last-commit/wanglab-cns/bidirectional-promoters?style=flat-square)
![GitHub issues](https://img.shields.io/github/issues/wanglab-cns/bidirectional-promoters?style=flat-square)
![GitHub pull requests](https://img.shields.io/github/issues-pr/wanglab-cns/bidirectional-promoters?style=flat-square)
![GitHub contributors](https://img.shields.io/github/contributors/wanglab-cns/bidirectional-promoters?style=flat-square)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/wanglab-cns/bidirectional-promoters?style=flat-square)

## Project Overview

A bidirectional promoter is a short stretch of DNA between two genes on opposite strands that start transcription
away from each other. One promoter controls two genes, so a single change there (for example, DNA methylation) can
affect both. One partner can be a protein-coding gene and the other a non-coding RNA such as a microRNA.

This project screens bidirectional promoters across common pediatric brain tumors and asks:

1. Which bidirectional promoters are altered in pediatric brain tumors, and is each partner gene silenced or overexpressed?
2. Are these changes specific to one tumor type, or shared across types?
3. Do they differ by molecular subtype, and are they linked to survival?

**Cohort:** Children's Brain Tumor Network (CBTN) tumors with whole-genome sequencing, RNA-seq, DNA methylation
and microRNA data. We start with six tumor types (low-grade glioma, ependymoma, medulloblastoma, ganglioglioma,
high-grade glioma and adamantinomatous craniopharyngioma) and then extend to adult datasets.

**Approach:**

- **Gene pairs:** build the list of bidirectional gene pairs from GENCODE v39, the version CBTN uses
- **Curation:** standardize clinical and molecular data into a `MultiAssayExperiment`, following the
  [ICB_Common](https://github.com/BHKLAB-DataProcessing/ICB_Common) curation standard and
  [ONCOTREE](https://oncotree.mskcc.org) tumor names
- **Screen:** test each promoter alteration against the expression of both partner genes, compared with
  single-direction promoters
- **Context:** compare results across tumor types and molecular subtypes, and with overall survival

## Set Up

### Prerequisites

Pixi is required to run this project.
If you haven't installed it yet, [follow these instructions](https://pixi.sh/latest/)

### Installation

```bash
# Clone the repository
git clone https://github.com/wanglab-cns/bidirectional-promoters.git
cd bidirectional-promoters

# Install dependencies via Pixi
pixi install
```

## Repository Structure

```
bidirectional-promoters/
├── config/              # Configuration files (dataset versions, cutoffs)
├── data/                # Not shared on GitHub; README in each folder
│   ├── rawdata/         # Downloaded data, never edited
│   ├── procdata/        # Curated and processed data
│   └── results/         # Analysis outputs
├── docs/                # Project documentation (MkDocs)
│   ├── README.md        # This page
│   ├── data_sources.md  # Every dataset: source, version, access date, citation
│   ├── usage.md         # How to set up and run the analysis
│   └── devnotes.md      # Dated design decisions
├── metadata/            # Dataset descriptions and mapping tables
├── workflow/
│   ├── notebooks/       # Exploration notebooks
│   └── scripts/         # Curation and analysis scripts, numbered in run order
└── pixi.toml            # Pixi environment
```

## Data

No data is stored in this repository, only code and documentation. Every dataset is listed in
[`docs/data_sources.md`](data_sources.md) with its source, version, access date and citation. Processed data that
can be shared publicly will be deposited on Zenodo, and the link will be added here.

## Documentation

Click [here](https://wanglab-cns.github.io/bidirectional-promoters) to view the full documentation.
