# Normal brain references for the CBTN tumor comparison

Checked 20 Oct 2026. I took every count, age and platform below from a GEO record, portal file or paper I opened. For GEO I parsed the per-sample metadata myself. Anything I could not confirm says "unverified". Every paper below was checked in PubMed for retractions and none are retracted. The one exception is Luo 2022 (Nature, PMID 36450980), which was retracted in 2025 and is left out.

## Recommendation

1. **Methylation.** No single dataset matches our platform, region and age, so build a pooled normal.
   - Process the raw IDATs with the same minfi pipeline we use for the tumors.
   - Keep only probes shared by 450K and EPIC v1.
   - Capper's controls (GSE90496) give the regions we can't get elsewhere: pons, hypothalamus and pituitary.
   - The pediatric sets give the age match:
     - Cerebellum: GSE278285 (24 children, 450K) and GSE299377 (EPIC).
     - Mixed regions: GSE215240 (13 pediatric EPIC normals, one each for pons, hypothalamus and pituitary).
     - Cortex: GSE74193 (76 children), GSE289191 (EPIC, 12 children) and PsychENCODE (450K, 13 children, 16 regions incl. temporal cortex and cerebellum).
2. **CBTN's own 23 EPIC normals.** Use them as a batch check, not as the main normal.
   - 21 of the 23 donors are tumor patients, mostly DMG.
   - No sample has a region label.
   - One sample looks like SHH tumor.
3. **Heidelberg v12.8 controls (GHGA, needs a data transfer agreement).** Apply only if optic-pathway LGG matters. It is the only source of normal optic nerve (CTRL_OPTIC).
4. **RNA.** Use GTEx v10 as the main normal.
   - It uses GENCODE v39, the same as our tumors.
   - It has cerebellum, hypothalamus, pituitary, spinal cord and cortex, all adult.
   - For age-matched layers add three sources:
     - PsychENCODE counts.
     - E-MTAB-6814 (open FASTQ, so it can be rerun with RSEM/v39).
     - GSE244124, the only pediatric pons and pituitary.
   - Put library type (poly-A vs rRNA-depleted) in the model.
5. **miRNA.** No HTG EdgeSeq normal brain exists in any public repository.
   - Fallback: GTEx v10 small RNA-seq (adult, includes pituitary and hypothalamus) plus PsychENCODE small RNA-seq (pediatric, 16 regions).
   - Use them for direction and rank only, not fold change.
   - The real fix is to ask CBTN/D3b to run HTG on a few of the normal-brain RNAs they already hold.

## Comparison table: methylation

| Dataset (first author, year) | Accession | Platform | Normal samples (n) | Brain regions (n per region) | Age range / pediatric n | Raw data? | Access | Best for | Caveats |
|---|---|---|---|---|---|---|---|---|---|
| CBTN/PBTA normals (OpenPedCan v15, 2025) | Kids First / OpenPedCan histologies v15 | EPIC v1 | 23 solid tissue + 1 blood | Not labelled. v12 classifier calls: CTRL_HEMI 6, CTRL_CBM 2, CTRL_CORPCAL 1, MB_SHH_2 1, no call 13 | 55 d–29.9 y; 20 under 18, only 1 under 5 y | IDATs held by CBTN (raw-data route unverified); betas open | Kids First / CBTN | Same-pipeline check; hemispheric HGG, ST-EPN | 21/23 donors are tumor patients (16 diffuse glioma, mostly DMG H3K27; 3 EPN; 1 GNT). The 55-d sample calls MB_SHH_2 (0.99) but the patient has EPN. Ask CBTN for regions. |
| Capper 2018 | GSE90496 | 450K | 72 normal-region (+24 inflammatory, 23 reactive) | PONS 12, CEBM 8, HYPTHAL 9, ADENOPIT 9, HEMI 13, WM 9, PINEAL 12 | Ages not on GEO (unverified); Cotter 2025 says the cerebellar controls are adult | IDATs yes | Open | Pons, hypothalamus, pituitary, all regions | PONS, HYPTHAL, ADENOPIT, PINEAL and WM are all FFPE. CEBM is 7 frozen + 1 FFPE; HEMI is 12 frozen + 1 FFPE. |
| Sill 2025/2026 (v12.8 reference) | GHGA GHGAS89861553411214 | 450K + EPIC + EPICv2 | n per class in Table S1 (unverified) | Classes: CTRL_ADENOPIT, CBM, CORPCAL, HEMI, HYPOTHAL, **OPTIC**, PIN, PONS, REACTIVE, BLOOD | Unverified | Raw, controlled | DTA, non-commercial use | Optic-pathway LGG (only optic control anywhere) | Slow to access; likely overlaps Capper. |
| Sturm 2023 | GSE215240 | EPIC | 13 | Hypothalamus 1, pituitary 1, pineal 1, pons 1, cerebellum+pons+medulla 1, cerebellum 1, 4th ventricle 1, temporal/parietal cortex 6 | 2–14 y; all 13 pediatric | IDATs yes | Open | ACP, DMG, PF-EPN, GG (pediatric, same platform) | FFPE surgical tissue from suspected-tumor cases, so it may be reactive or next to tumor. Only 1 sample per region. |
| Cotter 2025 | GSE299377 | EPIC v1 (4) + EPIC v2 (6) | 6 donors, 10 arrays | Cerebellar hemisphere cortex 6 | 0 mo (×2, premature to 10 d), 4 mo, 5 y, 16 y, 17 y | IDATs yes | Open | MB, cerebellar PA, PF-EPN | FFPE autopsy. Infant cerebellum looks like SHH-MB. Only donors 1–4 have EPIC v1. |
| Anne 2022 (GEO 2025) | GSE278285 | 450K | 45 controls (+30 autism) | Cerebellum 45 | 0–56 y; **24 under 18** (1 under 1 y, 11 aged 1–4, 6 aged 5–9, 6 aged 10–17) | IDATs yes | Open | MB, cerebellar PA, PF-EPN | Largest pediatric normal cerebellum set. Brain bank source not on GEO (unverified). Ages in whole years. |
| Li 2018 (PsychENCODE/BrainSpan) | development.psychencode.org; raw dbGaP phs000755 | 450K | 269 samples, 21 donors | 16 regions: ITC 19, A1C 18, OFC 18, VFC 18, STC 17, AMY 17, STR 17, MD thalamus 17, HIP 16, CBC 15, other cortex | 4 mo–40 y; 13 donors under 18 (159 samples); 8 pediatric CBC | Betas open (3.8 GB); IDATs controlled | Open / dbGaP | GG (temporal cortex), hemispheric HGG/LGG, thalamic HGG, MB | Healthy donors. No pons, hypothalamus or pituitary. |
| Jaffe 2016 | GSE74193 | 450K | 339 control donors on GEO (paper says 335); 450 arrays | DLPFC | Prenatal (36) to 85 y; **76 postnatal under 18** | IDATs yes | Open | Hemispheric HGG, ST-EPN, GG (cortex proxy) | Frontal cortex only. Use the bestQC flag. |
| Franklin 2025 | GSE289191 (bulk), GSE289184 (SATB2± sorted) | EPIC v1 | 107 bulk (95 fetal + 12 child); 87 sorted | Prefrontal cortex | Fetal 6–33 pcw; children 0–8 y (n=12) | IDATs yes | Open | Young hemispheric tumors | Frontal only; few children. |
| Nishitani 2025 | GSE239518 | EPIC | 11 typically developing (+13 maltreated) | Prefrontal cortex | 0–17.5 y; all pediatric | IDATs yes | Open | Hemispheric tumors | Forensic autopsies. The paper analysed 7 controls. |
| Hovestadt 2014 | GSE54880 | 450K | 8 | Cerebellum: 4 fetal, 4 adult | Fetal; adults 26–87 y; 0 children | **No IDATs** (beta + raw signal txt) | Open | MB backup | No children. |
| Smith 2022 | GSE209913 | EPIC | 14 | Fetal rhombic lip 7, EGL 7 (microdissected) | 14–21 pcw (paper) | IDATs yes | Open | MB Group 3/4 cell of origin | Progenitor zones, not mature tissue. |
| Kober 2018 | GSE115783 | 450K | 6 | Whole pituitary (autopsy) | Not on GEO | IDATs yes | Open | ACP | Adult pituitary-adenoma study. |
| Price 2016 | GSE69502 | 450K | 9 | Spinal cord | Fetal 14.5–20.4 wk | IDATs yes | Open | Spinal EPN (only option) | Fetal only. |
| Kobow 2020 / Jabari 2022 | GSE156374 | EPIC | 13 arrays (6 neocortex, 7 white matter) | Cortex, WM | Adult (mean ~30 y per papers) | IDATs yes | Open | GG (published comparator) | Adult FFPE autopsy. |
| Loyfer 2023 | GSE186458; raw EGA EGAS00001006791 | WGBS (sorted cells) | 14 brain | Cortex neurons 8, cerebellar neuron 1, neuron 1, oligodendrocytes 4 | 19–65 y | FASTQ controlled; processed open | Open / EGA | Cell-type deconvolution and markers | Adult, sorted cells, sequencing data. Not a bulk normal. |

## Comparison table: RNA

| Dataset (first author, year) | Accession | Platform | Normal samples (n) | Brain regions (n per region) | Age range / pediatric n | Raw data? | Access | Best for | Caveats |
|---|---|---|---|---|---|---|---|---|---|
| GTEx v10 (2024–25) | gtexportal.org; gs://adult-gtex/bulk-gex/v10; raw dbGaP phs000424 | Poly-A, unstranded, 76-bp PE; **GENCODE v39**; RNA-SeQC gene counts + RSEM transcripts | 3,234 brain (13 regions) + 313 pituitary | Cerebellum 266, cerebellar hemisphere 277, hypothalamus 257, pituitary 313, spinal cord C1 204, cortex 270, BA9 269, hippocampus 255, + 6 more | 20–79 y; 0 pediatric | FASTQ controlled | Counts/TPM open | All six (adult anchor) | Adult postmortem, poly-A only, no pons. |
| OpenPedCan v15 (Geng 2025) | AWS open-targets v15 / CAVATICA | GTEx v8 + PBTA, harmonized | 3 PBTA normal brain + GTEx v8 | PBTA normals have no region | PBTA: 218 d, 1,316 d, 4,968 d | Via Kids First | Open | Drop-in with our files | GTEx v8 here is not realigned. |
| Li 2018 (PsychENCODE/BrainSpan) | development.psychencode.org; raw dbGaP phs000755 | Poly-A, non-stranded, 76-bp SE; GENCODE v21 | 607 samples, 41 brains | 16 regions; CBC 33 | 8 pcw–40 y; 15 postnatal donors under 18 (222 samples); 14 pediatric CBC | Raw controlled | Counts open | MB, PF tumors, GG, hemispheric | No pons, hypothalamus or pituitary. Same donors as BrainSpan, so don't combine the two. |
| Cardoso-Moreira 2019 | E-MTAB-6814 | TruSeq stranded mRNA, 100-nt SE | 59 cerebellum/hindbrain + 55 brain | Cerebellum 59; cortex (DLPFC postnatal) 55 | Cerebellum: 33 prenatal, **17 under 18**, 9 adult; cortex 14 under 18 | FASTQ open | Open | MB, PF tumors (can rerun RSEM/v39) | 1–2 samples per age; early samples are whole hindbrain. |
| Griesinger 2024 (Donson lab) | GSE244124; SRA PRJNA1021260 | Poly-A, 150-bp PE | 16 | Cerebellum 6, pons 3, pituitary 2, frontal 2, parietal 1, temporal 1, cortex 1 | 0.2–9 y; all pediatric | FASTQ open | Open | **DMG (pons), ACP (pituitary)**, MB | 15/16 are autopsies from children with a brain-tumor history; 1 epilepsy resection. |
| Aldinger 2021 | dbGaP phs001908 | TruSeq RNA Access (capture) | 57 | Fetal whole cerebellum 13, RL 9, EGL 17, Purkinje layer 18 | 9–21 pcw | Controlled | dbGaP | MB cell of origin | Fetal; capture library. |
| HDBR (2016) | E-MTAB-4840 | Stranded PE | 632 | Cerebellum 36, spinal cord 42, diencephalon 42, temporal lobe 40, pons 2 | CS13–20 pcw | FASTQ open | Open | Fetal spinal cord, cerebellum | Fetal only; poly-A vs rRNA unverified. |
| Collado-Torres 2019 (BrainSeq 2) | eqtl.brainseq.org/phase2 | Ribo-Zero total RNA; GENCODE v25 | 300 DLPFC + 314 HIP controls | DLPFC, hippocampus | 14 pcw–96 y; pediatric n unverified | Raw access unverified | Processed open | Check for our rRNA-depleted libraries | Cortex and hippocampus only. |

## Comparison table: miRNA

**No HTG EdgeSeq normal brain exists in any public repository.** GEO has 53 HTG EdgeSeq series and none is normal brain tissue. OpenPedCan's 262 HTG samples are all tumors. Everything below is sequencing-based, so compare direction and rank only.

| Dataset (first author, year) | Accession | Platform | Normal samples (n) | Brain regions (n per region) | Age range / pediatric n | Raw data? | Access | Best for | Caveats |
|---|---|---|---|---|---|---|---|---|---|
| GTEx v10 small RNA (2025–26 release) | gs://adult-gtex/bulk-gex/v10/small-RNA | Small RNA-seq with UMIs; 1,979 miRNAs | 16,760 samples overall | Cerebellum 269, cerebellar hemisphere 288, hypothalamus 237, pituitary 314, spinal cord 136, cortex 238, BA9 249, hippocampus 237 | 20–79 y; 0 pediatric | Raw location unverified | Open | All six (adult anchor) | Not HTG; RNAcentral IDs need mapping to miRBase; no paper yet (unverified). |
| Li 2018 (PsychENCODE) | development.psychencode.org (smRNA count xlsx) | Small RNA-seq; 2,588 mature miRNAs | 278 samples, 22 postnatal brains | 16 regions; CBC 19 | 4 mo–40 y; 14 donors under 18 (157 samples); 12 pediatric CBC | dbGaP lists it; portal says "being uploaded" | Open | MB, PF tumors, GG, hemispheric | No pons, hypothalamus or pituitary. |
| BrainSpan miRNA (Allen; Ziats 2014) | download.alleninstitute.org/brainspan/MicroRNA | Platform unverified; 1,861 miRNAs | 216 samples, 18 donors | 16 structures | 4 mo–23 y (Allen API; the paper says 4 mo–19 y) | No | Open | Backup to PsychENCODE | Probably the same donors as PsychENCODE. |
| Keller 2022 (miRNATissueAtlas2) | GSE163534 | Small RNA-seq | 39 brain biopsies, 6 bodies | Cerebellum 2, pituitary 3, hypothalamus 1, spinal cord 2, temporal lobe 3, medulla 1 | Adult (ages unverified) | Yes | Open | Cross-check only | Tiny n. |

## Per tumor type

- **Medulloblastoma.** Methylation: GSE278285 (24 pediatric cerebellum), GSE299377 and PsychENCODE CBC; add GSE209913 for cell of origin. Infant cerebellum looks SHH-like, so don't treat it as clean normal. RNA: E-MTAB-6814 pediatric cerebellum, PsychENCODE CBC and GTEx cerebellum.
- **Low-grade glioma.** Cerebellar PA: same as MB. Optic/hypothalamic PA: Capper HYPTHAL plus the GSE215240 hypothalamus sample; the only optic control is Heidelberg CTRL_OPTIC (GHGA). Hemispheric: the cortex sets. RNA: GTEx cerebellum/hypothalamus plus PsychENCODE.
- **Ependymoma.** PF: Capper CEBM/PONS, GSE278285, and GSE215240 (4th ventricle, pons). ST: the cortex sets (GSE74193, GSE289191, PsychENCODE). Spinal: only fetal GSE69502. RNA: GTEx cerebellum/spinal cord plus GSE244124.
- **Ganglioglioma.** Methylation: PsychENCODE temporal cortex (ITC, STC and A1C: 31 pediatric samples) plus GSE215240 temporal cortex. GSE156374 is the adult comparator published GG papers used. RNA: PsychENCODE temporal cortex and GTEx cortex.
- **High-grade glioma.** Hemispheric: GSE74193, GSE289191, GSE239518, PsychENCODE and Capper HEMI/WM. Pons/DMG: Capper PONS (12) plus 2 GSE215240 pons samples; thalamic: PsychENCODE MD. RNA: GSE244124 pons (3), PsychENCODE and GTEx.
- **ACP.** Methylation: Capper ADENOPIT (9) + HYPTHAL (9), GSE115783 pituitary (6), and the GSE215240 pituitary + hypothalamus samples (1 each). RNA: GTEx pituitary/hypothalamus plus GSE244124 pituitary (2). No ACP methylation paper used a normal. ACP comes from Rathke's pouch epithelium, so pituitary and hypothalamus give local context, not the cell of origin.

## Gaps

- No HTG EdgeSeq normal brain, pediatric or adult.
- No public normal optic nerve methylation. CTRL_OPTIC is only in the Heidelberg set, behind a DTA.
- No postnatal spinal cord methylation on 450K or EPIC; only fetal (GSE69502). Spinal cord RNA and miRNA come only from adults (GTEx).
- Pediatric pons, hypothalamus and pituitary have only 1–3 samples each (GSE215240 EPIC; GSE244124 RNA). There is no pediatric pons, hypothalamus or pituitary miRNA.
- No healthy pediatric multi-region EPIC atlas. PsychENCODE is the closest, but it is 450K and has no pons or hypothalamus.
- Capper control ages are not public.
- Most pediatric normals are autopsy tissue from tumor patients (CBTN, GSE244124) or FFPE surgical tissue (GSE215240, Capper).
- Two datasets I could not get: Saratsis 2014 (7 pediatric brainstem, 450K) was never deposited. Mosella 2021's in-house pituitaries have a placeholder accession. Ask the authors for both.

## Practical notes

- Mixing 450K and EPIC: keep the shared probes, process everything from IDATs with one pipeline, and include platform and FFPE/frozen as covariates. Do not mix beta files from different pipelines.
- Region matters more than age for most promoters, and both matter in infants. Check neuron/glia proportion (e.g. with sorted references such as Loyfer 2023 or GSE289184) before calling a promoter hyper- or hypomethylated.
- RNA: GTEx v10 and our tumors share GENCODE v39. GTEx gene counts are from RNA-SeQC, not RSEM, so compare at TPM level or requantify.

## Also checked, not recommended

- GSE36278, GSE44684 and GSE77241: small 450K cerebellum/brain sets with no IDATs, likely the same DKFZ donors.
- GSE152561 (Lee 2024): 3 oxBS normals.
- GSE15745: 27K pons, too old.
- GSE53162: 6 control cerebellum, mostly adults.
- GSE131706: subventricular zone controls, 5 aged under 10.
- GSE73747, GSE58885 and GSE157908: fetal cortex only.
- GSE283928 (2025 PitNET): tumors only.
- Xena Toil TCGA-TARGET-GTEx: RSEM on GENCODE v23; a smaller GTEx subset.

## Sources

- Capper GSE90496: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90496
- Sill v12.8 (Cancer Cell 2026, doi:10.1016/j.ccell.2025.11.002); preprint medRxiv 10.1101/2025.05.28.25328344; full-text PDF https://pub.dzne.de/record/285055/files/DZNE-2026-00179.pdf ; GHGA https://data.ghga.de/browse?q=GHGAS89861553411214
- Sturm 2023 GSE215240: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215240 (PMID 36928815)
- Cotter 2025 GSE299377: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299377 (doi:10.1007/s00401-025-02907-1)
- Anne 2022 GSE278285: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278285 (doi:10.2217/epi-2022-0184)
- Li 2018 PsychENCODE: http://development.psychencode.org/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC6413317/ (doi:10.1126/science.aat7615)
- Jaffe 2016 GSE74193: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74193 (PMID 26619358)
- Franklin 2025 GSE289191 / GSE289184: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289191 (doi:10.1016/j.xgen.2025.101010)
- Nishitani 2025 GSE239518: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239518 (doi:10.1038/s41380-025-03236-1)
- Hovestadt 2014 GSE54880; Smith 2022 GSE209913; Kober 2018 GSE115783; Price 2016 GSE69502; Jabari 2022 GSE156374 (all at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<ID>)
- Loyfer 2023 GSE186458: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186458 (doi:10.1038/s41586-022-05580-6)
- OpenPedCan v15 histologies: https://s3.amazonaws.com/d3b-openaccess-us-east-1-prd-pbta/open-targets/v15/histologies.tsv (Geng 2025, PMC12402770)
- GTEx v10: https://gtexportal.org ; https://storage.googleapis.com/adult-gtex/bulk-gex/v10/small-RNA/description.tsv ; reference GTF gencode.v39 in gs://adult-gtex/references/v10/
- Cardoso-Moreira 2019 E-MTAB-6814: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6814 (PMC6658352)
- Griesinger 2024 GSE244124: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244124 (PMID 38334950)
- Aldinger 2021 (PMC8338761), HDBR E-MTAB-4840, BrainSeq Phase 2 (PMID 31174959)
- BrainSpan miRNA: https://download.alleninstitute.org/brainspan/MicroRNA/ ; miRNATissueAtlas2 GSE163534 (PMID 34570238)
- Luo 2022 retraction: PMID 36450980 (retraction note PMID 40634621)
