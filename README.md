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

## Curated analysis layout

The public-facing code is organised into five analytical blocks:

```text
analysis/
├── 01_bal_covid_validation/
├── 02_bcg_trained_immunity/
├── 03_cross_species_trm/
├── 04_luad_chemokine_correlations/
└── 05_pan_cancer_survival/
```

See [`analysis/README.md`](analysis/README.md) for the role and reproducibility status of each block, and [`docs/figure-map.md`](docs/figure-map.md) for the mapping to published panels.

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

## Relationship to the official code archive

The paper's code-availability statement points to the official scRNA-seq/TCR-seq workflow deposited by Victor Gourain:

- GitLab: https://gitlab.univ-nantes.fr/gourain-v-1/cxcr6_lungsepsis
- Zenodo: https://doi.org/10.5281/zenodo.10715057

That archive primarily covers the mouse single-cell RNA-seq and TCR-seq processing and downstream analysis. The present repository documents a different analytical layer: **human validation, cross-dataset integration and cancer-cohort analyses**.

## Data policy

**No study data are stored in this repository.**

The analyses rely on public or controlled-access resources documented in [`docs/data-sources.md`](docs/data-sources.md) and [`data/README.md`](data/README.md), including GEO, TCGA/UCSC Xena, EGA and the paper's source-data files.

Patient-level data, processed RDS objects, TCGA matrices, publisher source-data spreadsheets and controlled-access datasets are intentionally not duplicated here.

## Reproducibility status

These scripts originated within an active collaborative analysis environment in 2024. Some contain absolute local paths, references to objects generated in earlier preprocessing steps, legacy package interfaces or dependencies on data that must be obtained separately.

Accordingly, this repository should be read as **transparent scientific analysis code**, not a packaged one-command workflow. The curated code preserves the original scientific logic and filenames while making the published analysis path easier to navigate.

The complete original working layout is preserved on:

`archive/original-2024`

A pre-curation snapshot of the review branch is also preserved separately.

## Analytical narrative

The computational contribution represented here follows the study from mechanistic mouse findings toward human relevance:

**mouse mechanistic findings → public human datasets → cross-species tissue-resident T-cell signatures → cancer-cohort validation → survival relevance**

## Citation

If you use analyses from this repository, please cite the associated publication:

> Broquet A, Gourain V, Goronflot T, et al. Sepsis-trained macrophages promote antitumoral tissue-resident T cells. *Nature Immunology*. 2024;25:802–819. https://doi.org/10.1038/s41590-024-01819-8

For the official scRNA-seq/TCR-seq code release, also refer to Victor Gourain's Zenodo archive: https://doi.org/10.5281/zenodo.10715057

## Contact

**Debajyoti Sinha**  
Computational biology · single-cell transcriptomics · computational immunology
