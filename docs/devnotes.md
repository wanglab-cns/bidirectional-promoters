# Developer Notes

## Purpose of This Section

This section is for documenting technical decisions, challenges, and solutions encountered during your project. These notes are valuable for:

- Future you (who will forget why certain decisions were made)
- Collaborators who join the project later
- People coming from your publication who want to reproduce your work
- Anyone who might want to extend your research

## What to Document

### Design Decisions

Document important decisions about your project's architecture, algorithms, or methodologies:

``` markdown
## Choice of RNA-Seq Analysis Pipeline

[2025-04-25] We chose the kallisto over STAR pipeline for the following reasons:
    1. The CCLE dataset is very large, and kallisto is faster for quantifying large datasets
    2. GDSC used kallisto, so we can compare our results with theirs
```

### Technical Challenges

Record significant problems you encountered and how you solved them

``` markdown
## Sample Name Format Issue

[2025-04-25] We encountered a problem with sample name formats between the CCLE and GDSC datasets.
    The CCLE dataset uses "BRCA-XX-XXXX" format, while the GDSC dataset uses "BRCA-XX-XXXX-XX".
    We had to write a script to remove the last two characters from the sample names in the GDSC dataset.
```

### Dependencies and Environment

Document specific version requirements or compatibility issues:

``` markdown
## Critical Version Dependencies

[2025-04-25] SimpleITK 2.4.1 introduced a bug that flips images, so we froze version 2.4.0
```

## Best Practices

- Date your entries when appropriate
- Link to relevant code files or external resources
- Include small code snippets when helpful
- Note alternatives you considered and why they were rejected
- Document failed approaches to prevent others from repeating mistakes
- Update notes when major changes are made to the approach

---

# Project Decisions

## Scope: pan-brain, pediatric first

[2026-10-13] The screen covers common pediatric brain tumor types in CBTN rather than all cancer types. Pediatric
data comes first. Adult datasets (TCGA, CPTAC, GEO) are added once the pipeline works.

## GENCODE v39 for the gene pair list

[2026-10-13] CBTN data is built on GENCODE v39. An earlier pair list used GENCODE v50, but mixing versions causes
gene ID and coordinate mismatches downstream. The pair list is rebuilt on v39, and we count how many pairs are lost
compared with v50.

## Gene types

[2026-10-13] All gene types are kept, not only protein-coding genes, because microRNA and other non-coding partners
are part of the question.

## Clinical curation rules

[2026-10-13] We follow the [ICB_Common](https://github.com/BHKLAB-DataProcessing/ICB_Common) curation standard,
using [ICB_Van_Allen](https://github.com/BHKLAB-DataProcessing/ICB_Van_Allen-snakemake) as the example:

- Never delete columns. Key variables go first in the standard order, and every original column follows, unchanged.
- [2026-10-20] Column order (Farnoosh): the columns we care about go first on the left, then the rest. In `CLIN.csv` that is the lab's key columns, then the columns we use (tumor type, subtype, ONCOTREE, same-sample flag, data types, demographics, location, resection, EFS), then every other original column. NB01's tables follow the same rule.
- No imputation. Missing values stay `NA`. Variables we don't have (drug, response, PFS, stage, histology) are `NA`.
- Survival time is in months; the event is 1 = death, 0 = alive. Overall survival only for now. Event-free
  survival waits until its event definition is clear, but its original columns are kept.
- Sex is `F`/`M`. Age at diagnosis is `age`.
- Tumor type and subtype names follow ONCOTREE. All six main tumor types sit at ONCOTREE level 3.
- Each data release is versioned in its name, for example `CBTN_2026`.

## Analysis cohort

[2026-10-20] The main set for analyses that combine data types is the 186 patients with WGS, RNA-seq, methylation and
miRNA from the same tumor sample (`core_mirna` in NB01). All 203 miRNA patients stay in the curated tables, so analyses
that need fewer data types (for example survival) can use the other 17. Those 17 hold 7 of the 36 deaths.

## Repository template, kept in full

[2026-10-20] Farnoosh said on 13 Oct that the full BHK template is more than this project needs (data, docs and workflow
folders are enough). We keep the full template by choice; it has those folders, and it is the one she used for
PredictIO-UV-Dist.

## No data on GitHub

[2026-10-13] This repository is public and holds code and documentation only. Publicly shareable processed data
will go to Zenodo, with the link in the README.

## Repository template

[2026-10-20] Set up from [bhklab-project-template](https://github.com/bhklab/bhklab-project-template) v0.14.1 with
copier. The template's setup tasks were skipped because the repository already existed, so `pixi.lock` is created on
the first `pixi install`.
