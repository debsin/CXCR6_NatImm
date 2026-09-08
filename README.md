# CXCR6_NatImm

[![DOI](https://img.shields.io/badge/Nature%20Immunology-10.1038%2Fs41590--024--01819--8-blue)](https://doi.org/10.1038/s41590-024-01819-8)

**Computational analyses supporting the human and cross-dataset validation components of a collaborative Nature Immunology study on sepsis-trained macrophages, CXCR6+ tissue-resident T cells and cancer immune surveillance.**

This repository contains analysis scripts developed in support of the human and cross-dataset validation components of:

**Broquet A, Gourain V, Goronflot T, et al.**  
*Sepsis-trained macrophages promote antitumoral tissue-resident T cells.*  
**Nature Immunology. 2024;25:802–819.**  
DOI: [10.1038/s41590-024-01819-8](https://doi.org/10.1038/s41590-024-01819-8)

The analyses here principally correspond to **Figure 8** and **Extended Data Figure 10** in the published article. They focus on human validation, cross-species projection, public-dataset integration, cancer-cohort analysis and survival relevance.

This repository is **complementary** to the study's official scRNA-seq and TCR-seq code archive maintained by Victor Gourain. It does not represent the complete analytical workflow of the collaborative study.

## Contribution scope

The computational work represented here includes:

- analysis of publicly available bronchoalveolar lavage single-cell datasets from healthy, moderate-infection and severe/sepsis COVID-19 cohorts;
- evaluation of CXCR6 expression in T-cell and tissue-resident T-cell compartments;
- analysis of CXCL16, CCL2 and CCL7 expression in macrophage compartments;
- trained-immunity signature analysis following BCG vaccination;
- definition and cross-species projection of human and murine tissue-resident T-cell signatures;
- comparison of CXCR6 expression with tissue-resident T-cell signatures in human datasets;
- TCGA cancer analyses linking TRM, macrophage and chemokine signatures;
- CXCR6–CCR2 coexpression analysis and survival stratification across cancer cohorts.

For a panel-level map, see [`docs/figure-map.md`](docs/figure-map.md).

## Relationship to the official code archive

The paper's code-availability statement points to the official scRNA-seq/TCR-seq workflow deposited by Victor Gourain:

- GitLab: https://gitlab.univ-nantes.fr/gourain-v-1/cxcr6_lungsepsis
- Zenodo: https://doi.org/10.5281/zenodo.10715057

That archive primarily covers the mouse single-cell RNA-seq and TCR-seq processing and downstream analysis. The present repository documents a different analytical layer: **human validation, cross-dataset integration and cancer-cohort analyses**.

## Data policy

No patient-level data, processed RDS objects, TCGA matrices, publisher source-data spreadsheets or controlled-access datasets should be stored in this repository.

The analyses instead rely on public or controlled-access sources documented in [`docs/data-sources.md`](docs/data-sources.md), including GEO, TCGA/UCSC Xena, EGA and the paper's own source-data files.

## Reproducibility status

The original 2024 scripts were written as analysis code within an active collaborative project. Several contain absolute local paths, references to objects created in earlier sessions, legacy Seurat slot access or dependencies on datasets that must be obtained separately.

The curated public version should therefore distinguish between:

- analyses reproducible directly from public data;
- analyses reproducible after obtaining controlled-access data;
- documented analysis code preserved for transparency but not yet packaged as a one-command workflow.

The unmodified original repository state is preserved on the branch:

`archive/original-2024`

This branch, `public-review`, is used to prepare a clearer public-facing version without altering the original history.

## Repository direction

The goal of the curated version is not to claim ownership of the full study, but to make a specific computational contribution visible and traceable:

**mouse mechanistic findings → public human datasets → cross-species tissue-resident T-cell signatures → cancer-cohort validation → survival relevance**

## Citation

If you use analyses from this repository, please cite the associated publication:

> Broquet A, Gourain V, Goronflot T, et al. Sepsis-trained macrophages promote antitumoral tissue-resident T cells. *Nature Immunology*. 2024;25:802–819. https://doi.org/10.1038/s41590-024-01819-8

For the official scRNA-seq/TCR-seq code release, also refer to Victor Gourain's Zenodo archive: https://doi.org/10.5281/zenodo.10715057

## Contact

**Debajyoti Sinha**  
Computational biology · single-cell transcriptomics · computational immunology
