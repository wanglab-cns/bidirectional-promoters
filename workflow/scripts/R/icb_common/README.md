# ICB_Common (BHK lab)

Curation functions from the BHK lab, copied without changes from
[BHKLAB-DataProcessing/ICB_Common](https://github.com/BHKLAB-DataProcessing/ICB_Common/tree/308a68596c9c5aa00fd4278d5605c068c6198aad/code)
(the old name BHKLAB-Pachyderm/ICB_Common redirects there), commit `308a685` (4 Jun 2026), folder `code/`. The `code/`
files last changed on 1 Dec 2022.

Our notebooks `source()` these local copies instead of the GitHub URLs, so the code we ran stays fixed.
Any change we need is made in a later, separate commit with a comment that says what changed and why.

| File | Used in | What it does |
|---|---|---|
| `format_clin_data.R` | NB02 | Appends the original clinical columns to the curated ones |
| `annotate_tissue.R` | NB02 | Adds `tissueid` from `curation_tissue.csv`, survival unit and type, standard column order |
| `format_clin_cased_sequenced.R` | NB04 | Saves `CLIN.rds` and `cased_sequenced.rds` |
| `format_se.R`, `Create_SummarizedExp.R` | NB04 | Build SummarizedExperiments, rename clinical columns (e.g. `t.os` to `survival_time_os`) |
| `get_exp_se.R`, `get_snv_se.R`, `get_cna_se.R` | NB03, NB04 | Expression, mutation and copy-number SummarizedExperiments |
| `SNV_Functions.R`, `CNA_Functions.R`, `get_snv_feat.R`, `get_cna_feat.R` | NB04 | Mutation and copy-number features (TMB, CIN, CNA counts) |
| `get_multiassay.R`, `get_MultiAssayExp.R` | NB04 | Combine the SummarizedExperiments into a MultiAssayExperiment |

Known changes we will need (made later, not here):

- `get_*_se.R` and `get_MultiAssayExp.R` still `source()` code from GitHub URLs; NB04 will point them at these local copies.
- The SNV builder keeps coding variants only. We need promoter (non-coding) mutations too.
- `get_MultiAssayExp` assumes one sample per patient. Our patients can have several tumor samples, so it needs a sample map.
