# Metadata Directory

Small descriptive files that the curation needs, for example:

- the mapping from CBTN tumor labels to ONCOTREE names and codes, with the ONCOTREE version used
- column dictionaries for curated tables

No patient-level data goes here.

## Files

| File | Content |
|---|---|
| `curation_tissue.csv` | Tumor type to tissue (ONCOTREE level 1) and type code, read by the lab's `annotate_tissue()` (NB02) |
| `oncotree_mapping.csv` | Ordered rules from CBTN tumor type, molecular subtype and pathology text to the most detailed ONCOTREE code (`oncotree_2025_10_03`); the first matching rule wins (NB02) |
