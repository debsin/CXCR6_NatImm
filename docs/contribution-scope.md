# Contribution scope

This repository documents a **specific computational contribution** to a broader collaborative study.

## What this repository represents

The code here supports the human and cross-dataset validation layer of the Nature Immunology study:

> Broquet A, Gourain V, Goronflot T, et al. *Sepsis-trained macrophages promote antitumoral tissue-resident T cells.* Nature Immunology. 2024;25:802–819. https://doi.org/10.1038/s41590-024-01819-8

The analyses represented here principally contributed to Figure 8 and Extended Data Figure 10 and include:

- human BAL single-cell validation;
- human and murine tissue-resident T-cell signature comparison;
- cross-species projection;
- trained-immunity signature analysis;
- TCGA chemokine/macrophage/TRM correlations;
- CXCR6–CCR2 coexpression analyses;
- cancer survival analyses across multiple cohorts.

## Scope boundary

This is not the complete computational workflow for the paper. The mouse scRNA-seq and TCR-seq processing and downstream-analysis workflow is archived separately at:

- https://gitlab.univ-nantes.fr/gourain-v-1/cxcr6_lungsepsis
- https://doi.org/10.5281/zenodo.10715057

The present repository focuses on the human validation, cross-dataset integration and cancer-cohort analyses.

## Attribution

The paper's published contribution statement does not separately enumerate individual computational contributions for the Figure 8 / Extended Data Figure 10 analyses.

For that reason, this repository uses the formulation:

> This repository contains analysis scripts developed in support of the human and cross-dataset validation components of the collaborative study.

It does not claim sole authorship of Figure 8 or of the complete computational workflow.
