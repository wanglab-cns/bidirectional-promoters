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
- **Version**: File manifests downloaded 2026-10-06, labelled `CBTN_2026` in this project. This release has more
  patients than OpenPedCan v15.
- **URL**: https://cavatica.sbgenomics.com
- **Access Method**: CAVATICA account and acceptance of the CBTN data access terms; file manifests exported from the
  CAVATICA projects
- **Access Date**: 2026-10-06
- **Data Format**: TSV manifests; IDAT (methylation)
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
- **Version**: `oncotree_latest_stable`. The exact version is recorded in `metadata/` when the mapping is built.
- **URL**: https://oncotree.mskcc.org (API: https://oncotree.mskcc.org/api/)
- **Access Method**: API
- **Citation**: Kundra R, ... Schultz N. OncoTree: A Cancer Classification System for Precision Oncology.
  *JCO Clinical Cancer Informatics* (2021). [doi:10.1200/CCI.20.00108](https://doi.org/10.1200/CCI.20.00108)

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
- **Input Data**: GENCODE v39
- **Processing Scripts**: `workflow/scripts/` (to be added)

### Curated clinical data

- **Name**: Curated clinical table for CBTN patients with microRNA data
- **Creation Method**: Clinical curation following the
  [ICB_Common](https://github.com/BHKLAB-DataProcessing/ICB_Common) standard, using
  [ICB_Van_Allen](https://github.com/BHKLAB-DataProcessing/ICB_Van_Allen-snakemake) as the worked example
- **Input Data**: CBTN / OpenPedCan v15, CBTN on CAVATICA, PedcBioPortal
- **Processing Scripts**: `workflow/scripts/` (to be added)

---

## 3. Data Dictionary

### Curated clinical table (draft)

Key variables come first, in this order. Every original column is kept after them, unchanged. Missing values
stay `NA`, and nothing is imputed or removed.

| Column Name | Data Type | Description | Units | Possible Values |
|---|---|---|---|---|
| patient | string | CBTN participant ID | N/A | `PT_XXXXXXXX` |
| sex | string | Reported sex | N/A | F, M, NA |
| age | numeric | Age at diagnosis | years | ≥ 0 |
| primary | string | Tumor type (ONCOTREE name) | N/A | e.g. Medulloblastoma |
| histo | string | Histology | N/A | NA (not available) |
| tissueid | string | Tissue (ONCOTREE level 1) | N/A | CNS/Brain |
| stage | string | Stage | N/A | NA (not used for these tumors) |
| response.other.info, recist, response | string | Treatment response | N/A | NA (no response data) |
| drug_type, treatmentid | string | Treatment | N/A | NA (no drug data) |
| dna | string | DNA data type | N/A | wgs, NA |
| rna | string | RNA data type | N/A | e.g. tpm, NA |
| t.pfs, pfs | numeric | Progression-free survival | months; 0/1 | NA (event-free survival kept in the original columns for now) |
| t.os | numeric | Overall survival time | months | ≥ 0 |
| os | integer | Overall survival event | N/A | 1 = death, 0 = alive |
| survival_unit | string | Unit of survival time | N/A | month |
| survival_type | string | Survival endpoints available | N/A | OS |

---

## Best Practices

- Store raw data in `data/rawdata/` and never modify it
- Store processed data in `data/procdata/`; all code used to generate it goes in `workflow/scripts/`
- Document all processing steps
- Track data provenance (where data came from and how it was modified)
- Respect data usage agreements and licenses. CBTN patient-level data must never be pushed to GitHub.
