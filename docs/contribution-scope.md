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

## What this repository does not represent

This is not the complete computational workflow for the paper.

The official mouse scRNA-seq and TCR-seq processing and downstream-analysis code was maintained separately by Victor Gourain and is archived at:

- https://gitlab.univ-nantes.fr/gourain-v-1/cxcr6_lungsepsis
- https://doi.org/10.5281/zenodo.10715057

The present repository is therefore complementary rather than redundant.

## Attribution

The paper's published contribution statement lists Debajyoti Sinha among the authors who performed experiments and contributed to interpretation of results and manuscript revision. It does not separately enumerate individual computational contributions for the Figure 8 / Extended Data Figure 10 analyses.

For that reason, the repository should use language such as:

> This repository contains analysis scripts developed in support of the human and cross-dataset validation components of the collaborative study.

It should avoid claims such as:

> Code for Figure 8 developed solely by Debajyoti Sinha.

unless a more specific attribution has been explicitly agreed by the collaborators.

## Professional positioning

The scientific value of this contribution lies in connecting mechanistic mouse findings to human evidence:

**mouse immunology → public human single-cell datasets → cross-species signatures → cancer cohorts → survival relevance**

This is appropriately described as computational immunology and cross-dataset translational validation.
