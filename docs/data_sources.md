# Data Sources

## Overview

This page lists every dataset used in the project: where it came from, which version, when it was downloaded, and
how to cite it. Raw files go in `data/rawdata/<source>/` and are never edited. Nothing in `data/` is pushed to GitHub.

---

## 1. External Data Sources

### CBTN / OpenPedCan v15 (open access)

- **Name**: Open Pediatric Cancer (OpenPedCan) project, release v15, which includes the Children's Brain Tumor
  Network (CBTN) / Pediatric Brain Tumor Atlas (PBTA) cohort
- **Version**: v15 (release date 2024-03-01)
- **URL**: https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/open-targets/v15/
- **Access Method**: Direct download over HTTPS, checked against the release `md5sum.txt` and `methyl-md5sum.txt`
- **Access Date**: 2026-10-06
- **Data Format**: TSV (`histologies.tsv`, independent-specimen lists, gene ID map, EPIC probe annotation),
  RDS (`mirna-expression-counts.rds`, expression matrices)
- **Citation**:
    - Geng Z, ... Rokita JL. The Open Pediatric Cancer Project. *GigaScience* (2025).
      [doi:10.1093/gigascience/giaf093](https://doi.org/10.1093/gigascience/giaf093)
    - Shapiro JA, ... OpenPBTA: The Open Pediatric Brain Tumor Atlas. *Cell Genomics* (2023).
      [doi:10.1016/j.xgen.2023.100340](https://doi.org/10.1016/j.xgen.2023.100340)
- **License**: See the OpenPedCan data release terms

### CBTN harmonized data on CAVATICA (latest release)

- **Name**: CBTN harmonized data and CBTN-Methylation source data (raw EPIC IDAT files)
- **Version**: Project `cavatica/pbta-cbttc` ("PBTA-CBTN Open Access"). File manifests first exported 2026-10-06;
  manifests and files for our cohort downloaded 2026-10-20 by `workflow/scripts/00_download_cbtn_2026.py`
  (1,888 files, 8.7 GB: RSEM gene results, consensus somatic MAF, copy-number segments, EPIC IDATs, plus the
  non-tumor methylation samples). Labelled `CBTN_2026` in this project. This release has more patients than
  OpenPedCan v15. Its RNA files use GENCODE v39.
- **URL**: https://cavatica.sbgenomics.com
- **Access Method**: CAVATICA account and acceptance of the CBTN data access terms; file manifests exported from the
  CAVATICA projects
- **Access Date**: 2026-10-06 (manifests), 2026-10-20 (files)
- **Data Format**: TSV manifests; RSEM results; MAF; SEG; IDAT (methylation)
- **Citation**: Chen F, ... Creighton CJ. The DNA methylome of pediatric brain tumors appears shaped by structural
  variation and predicts survival. *Nature Communications* (2024).
  [doi:10.1038/s41467-024-51276-y](https://doi.org/10.1038/s41467-024-51276-y)
- **License**: CBTN data access terms. Patient-level data must not be shared outside the approved users.

### PedcBioPortal: PBTA (Provisional)

- **Name**: Pediatric Brain Tumor Atlas, all samples (`pbta_all`)
- **Version**: Portal exports downloaded 2026-10-06
- **URL**: https://pedcbioportal.kidsfirstdrc.org
- **Access Method**: Sample and clinical tables exported from the portal
- **Access Date**: 2026-10-06
- **Data Format**: TSV
- **Citation**: as for CBTN / OpenPedCan above

### GENCODE v39 gene annotation

- **Name**: GENCODE human gene annotation, release 39 (GRCh38). CBTN data is built on this release, so we use it
  for the bidirectional pair list.
- **Version**: Release 39
- **URL**: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz
- **Access Method**: Direct download, checked against the release `MD5SUMS`
- **Access Date**: 2026-10-06
- **Data Format**: GTF
- **Citation**: Frankish A, ... Flicek P. GENCODE 2021. *Nucleic Acids Research* (2021).
  [doi:10.1093/nar/gkaa1087](https://doi.org/10.1093/nar/gkaa1087)

### ONCOTREE

- **Name**: OncoTree cancer type classification, used to name tumor types and subtypes (levels 1 to 3)
- **Version**: `oncotree_2025_10_03` (the latest stable version on 2026-10-20)
- **URL**: https://oncotree.mskcc.org (API: https://oncotree.mskcc.org/api/tumorTypes?version=oncotree_2025_10_03)
- **Access Method**: the BHK lab R package [AnnotationGx](https://github.com/bhklab/AnnotationGx) (`getOncotreeTumorTypes()`, version 0.99.3, GitHub commit `adb3838`, the last that runs on R 4.5), saved as `data/rawdata/oncotree_2025_10_03/tumor_types_annotationgx.csv`; the API JSON is kept next to it
- **Access Date**: 2026-10-20
- **Citation**: Kundra R, ... Schultz N. OncoTree: A Cancer Classification System for Precision Oncology.
  *JCO Clinical Cancer Informatics* (2021). [doi:10.1200/CCI.20.00108](https://doi.org/10.1200/CCI.20.00108)

### UCSC CpG islands (hg38)

- **Name**: UCSC `cpgIslandExt` track, GRCh38, used to mark pairs whose shared promoter overlaps a CpG island
- **URL**: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz (last modified 2022-10-19)
- **Access Date**: 2026-10-20, saved in `data/rawdata/ucsc_hg38_cpg/`

### GENCODE v50 gene annotation (comparison only)

- **Name**: GENCODE release 50, used only to count the pairs gained and lost against v39
- **URL**: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz
- **Access Date**: 2026-10-20 (MD5 checked), saved in `data/rawdata/gencode_v50/`

### Normal brain references (candidates)

Reviewed on 2026-10-20; none downloaded yet. See [normal_reference.md](normal_reference.md).

### Adult datasets (planned, second phase)

Only clinical and annotation files have been downloaded so far (2026-10-09); no omics data yet.

| Dataset | What we use it for | Source |
|---|---|---|
| TCGA-LGG, TCGA-GBM | Adult glioma: RNA-seq, miRNA-seq, methylation | [GDC](https://portal.gdc.cancer.gov), Data Release 46.0 |
| CPTAC-3 GBM | Adult glioblastoma with EPIC methylation | [GDC](https://portal.gdc.cancer.gov), Data Release 46.0 |
| TCGA Pan-Cancer Clinical Data Resource | Curated TCGA survival | Liu J et al., *Cell* (2018), [doi:10.1016/j.cell.2018.02.052](https://doi.org/10.1016/j.cell.2018.02.052) |
| Capper 2018 (GSE90496, GSE109379) | Methylation of all six types, with age | Capper D et al., *Nature* (2018), [doi:10.1038/nature26000](https://doi.org/10.1038/nature26000) |
| Cavalli 2017 (GSE85218) | Medulloblastoma expression and methylation, including adults | Cavalli FMG et al., *Cancer Cell* (2017), [doi:10.1016/j.ccell.2017.05.005](https://doi.org/10.1016/j.ccell.2017.05.005) |
| GSE243240, GSE64415 | Ependymoma | [GEO](https://www.ncbi.nlm.nih.gov/geo/) |
| GSE239695, GSE294056 | Adamantinomatous craniopharyngioma | [GEO](https://www.ncbi.nlm.nih.gov/geo/) |

---

## 2. Internal/Generated Data

Each entry is added here when the file is first made.

### Bidirectional gene pairs (GENCODE v39)

- **Name**: Genome-wide list of bidirectional gene pairs
- **Creation Method**: Head-to-head gene pairs on opposite strands, with transcription start sites within a set
  distance, built from the GENCODE v39 GTF
- **Input Data**: GENCODE v39, UCSC CpG islands
- **Processing Scripts**: `workflow/scripts/01_bidirectional_pairs.py` (run with `--canonical-only`)
- **Output**: `data/results/pairs_GENCODE_v39/` (not on GitHub)

### Curated clinical data

- **Name**: Curated clinical table for CBTN patients with microRNA data
- **Creation Method**: Clinical curation following the
  [ICB_Common](https://github.com/BHKLAB-DataProcessing/ICB_Common) standard, using
  [ICB_Van_Allen](https://github.com/BHKLAB-DataProcessing/ICB_Van_Allen-snakemake) as the worked example
- **Input Data**: OpenPedCan v15 `histologies.tsv` (survival, age, sex, subtype); newest CBTN manifests (check of sex);
  NB01 outputs (cohort, data types per patient); ONCOTREE
- **Processing**: `workflow/notebooks/02_clinical.ipynb`, using the lab functions in `workflow/scripts/R/icb_common/`
- **Output**: `data/procdata/02_clinical/CLIN.csv` (not on GitHub); `metadata/curation_tissue.csv`

---

## 3. Data Dictionary

### Curated clinical table (`CLIN.csv`)

One row per patient, `;`-separated, written as in the lab code. The lab's key columns come first, in the lab's order,
then our extra columns, then every original v15 column of the patient's first miRNA tumor. Missing values stay `NA`;
nothing is imputed. A `;` inside a value is written as ` | ` so the file stays readable. Rules: `devnotes.md`
(clinical curation rules, 13 Oct).

| Column Name | Data Type | Description | Units | Possible Values |
|---|---|---|---|---|
| patient | string | CBTN participant ID | N/A | `PT_XXXXXXXX` |
| sex | string | Reported sex (v15 `reported_gender`) | N/A | F, M |
| age | numeric | Age at diagnosis (earliest over the patient's miRNA samples) | years | 0 to 31 |
| primary | string | Tumor type of the first miRNA tumor, as its ONCOTREE name | N/A | e.g. Medulloblastoma, Low-Grade Glioma, NOS |
| histo | string | Histology | N/A | NA (not available) |
| tissueid | string | Tissue (ONCOTREE), from `metadata/curation_tissue.csv` | N/A | CNS/Brain |
| treatmentid | string | Treatment | N/A | empty (no drug data) |
| stage | string | Stage | N/A | NA (not used for these tumors) |
| response.other.info, recist, response | string | Treatment response | N/A | NA (no response data) |
| drug_type | string | Treatment type | N/A | NA |
| dna | string | DNA data type | N/A | wgs, NA |
| rna | string | RNA data type (as in the lab's 2026 ICB repos) | N/A | rnaseq, NA |
| t.pfs, pfs | numeric | Progression-free survival | months; 0/1 | NA (CBTN reports event-free survival, kept in the original columns `EFS_days`, `EFS_event_type`) |
| t.os | numeric | Overall survival time (`OS_days` / 30.5) | months | ≥ 0 |
| os | integer | Overall survival event | N/A | 1 = deceased, 0 = living |
| survival_unit | string | Unit of survival time | N/A | month |
| survival_type | string | Survival endpoints available | N/A | OS |
| tumor_type | string | CBTN tumor type(s) of the patient's miRNA samples (`cancer_group`) | N/A | the six types; one patient has two, joined with ` \| ` |
| rna_info, dna_info | string | What the RNA and DNA data are | N/A | RSEM TPM and expected counts; consensus somatic MAF, CNVkit segments; NA |
| methyl, methyl_info | string | Methylation data type and processing | N/A | epic; processing NA until decided |
| mirna, mirna_info | string | miRNA data type and processing | N/A | htg; raw counts |
| race, ethnicity | string | As reported in v15 ("Reported Unknown", "Not Available" kept as written) | N/A | |
| oncotree_code, oncotree_name, oncotree_level, oncotree_level1 ... oncotree_level5 | string | Most detailed ONCOTREE code the data supports (Farnoosh's rule), from the rules in `metadata/oncotree_mapping.csv`, with its name, level, and the name at each level (1 = tissue, 2 = family, 3 = type, 4-5 = subtypes; NA below the code's level). A patient whose samples map to different codes keeps all, joined with ` / ` | N/A | e.g. MBLG4, Medulloblastoma, Group 4, 5 |
| sex_newest_release | string | Sex in the newest CBTN release, as a check | N/A | F, M, NA |
| n_mirna_samples | integer | miRNA samples of the patient in the six types | N/A | 1 to 3 |
| (original columns) | | Every v15 `histologies.tsv` column of the patient's first miRNA tumor, except those turned into key columns | | |

---

## Best Practices

- Store raw data in `data/rawdata/` and never modify it
- Store processed data in `data/procdata/`; all code used to generate it goes in `workflow/scripts/`
- Document all processing steps
- Track data provenance (where data came from and how it was modified)
- Respect data usage agreements and licenses. CBTN patient-level data must never be pushed to GitHub.
