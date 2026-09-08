# Data sources

This repository should not contain patient-level datasets, processed RDS objects, TCGA matrices, controlled-access data or publisher source-data workbooks.

The analyses rely on external datasets and published resources that should be obtained from their original repositories.

## Study data and official code

- Nature Immunology article: https://doi.org/10.1038/s41590-024-01819-8
- Single-cell sequencing accession: ENA `PRJEB52332`
- Official scRNA-seq/TCR-seq code archive: https://doi.org/10.5281/zenodo.10715057
- Victor Gourain's GitLab repository: https://gitlab.univ-nantes.fr/gourain-v-1/cxcr6_lungsepsis

## Human BAL / COVID-19 single-cell data

- GEO: `GSE145926`

Used for analyses of CXCR6 expression in T-cell/TRM compartments and CXCL16, CCL2 and CCL7 expression in macrophages across healthy, moderate-infection and severe/sepsis groups.

## Trained-immunity / BCG datasets

- GEO: `GSE124220`
- Related study resources include `GSE124218`

Used for projection/evaluation of trained macrophage or trained-immunity signatures after BCG vaccination.

## Human tissue-resident T-cell reference datasets

The analyses use published TRM signatures and supplementary tables from prior human lung, tumor and dermis studies. These source tables should be retrieved from the relevant publications rather than committed here.

## Human lung-tumor T-cell data

- EGA: `EGAS00001004707`

This is a controlled-access resource. Analyses depending on it cannot be reproduced without appropriate data access.

## Cancer cohorts

TCGA RNA-seq and survival data were used for lung adenocarcinoma and additional tumor types.

Original analyses referenced UCSC Xena / GDC-derived TCGA matrices and associated survival information. Users should obtain current versions from the official data portals rather than from this repository.

- UCSC Xena: https://xenabrowser.net/datapages/
- NCI Genomic Data Commons: https://portal.gdc.cancer.gov/

## Publisher source data

The Nature Immunology article provides source-data workbooks for Figure 8 and Extended Data Figure 10. These are useful for checking the final published statistical values and panel composition but should not be duplicated in this repository.

## Reproducibility classification

Analyses in the curated repository should be labelled as one of:

1. **Public-data reproducible** — all required data can be downloaded openly.
2. **Controlled-access reproducible** — code is provided, but the user must obtain permission for one or more datasets such as EGA resources.
3. **Documented analysis code** — original collaborative analysis is preserved for transparency but is not yet packaged into a fully standalone workflow.
