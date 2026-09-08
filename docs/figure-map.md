# Figure map

This document maps the principal computational analyses in this repository to the published Nature Immunology article.

## Published Figure 8

| Panel | Published analysis | Curated repository material |
| --- | --- | --- |
| 8a | CXCR6 expression in BAL T-cell and tissue-resident T-cell compartments across healthy, moderate infection and severe/sepsis COVID-19 | `analysis/01_bal_covid_validation/NatImm_Figure9.R` |
| 8b | CXCL16, CCL7 and CCL2 expression in BAL macrophages | `analysis/01_bal_covid_validation/NatImm_Figure9.R` |
| 8c | Trained alveolar-macrophage signature after BCG vaccination | `analysis/02_bcg_trained_immunity/process_BCG.R` |
| 8d | Projection of a human lung-tumor TRM signature onto mouse lung scRNA-seq | `analysis/03_cross_species_trm/compare_sig.R`, `compare_sig_human.R`, `plots_TRM.R` |
| 8e | Enrichment of the murine CXCR6/TRM signature in human tumor TRM versus non-TRM cells | `analysis/03_cross_species_trm/compare_sig.R`, `compare_sig_human.R`, `plots_TRM.R` |
| 8f | LUAD correlations between TRM signature, CXCL16, macrophage markers and CCR2 ligands | `analysis/04_luad_chemokine_correlations/panel_TCGA.R` |
| 8g | CXCR6–CCR2 coexpression and survival in LUAD | `analysis/05_pan_cancer_survival/NatImm_Figure9_G-J.R` |
| 8h–j | CXCR6–CCR2 coexpression and survival across additional TCGA cancers | `analysis/05_pan_cancer_survival/NatImm_Figure9_G-J.R` |

## Extended Data Figure 10

| Panel | Published analysis | Curated repository material |
| --- | --- | --- |
| 10a | CXCL16, CCL7 and CCL2 expression in BAL macrophages | `analysis/01_bal_covid_validation/NatImm_Figure9.R` |
| 10b–c | Human tumor immunostaining | Experimental imaging; not a computational analysis represented here |
| 10d | Correlation of CXCR6 with human lung-tumor, lung and dermis TRM signatures | `analysis/04_luad_chemokine_correlations/panel_TCGA.R` |
| 10e | Projection of human dermis TRM signature onto mouse scRNA-seq | `analysis/03_cross_species_trm/compare_sig.R`, `compare_sig_human.R`, `plots_TRM.R` |
| 10f | Definition of a murine CXCR6 tissue-resident T-cell signature | `analysis/03_cross_species_trm/compare_sig.R`, `plots_TRM.R` |

## Naming note

Several preserved filenames still refer to `Figure9`. During manuscript development, the human-validation figure was numbered differently; it became **Figure 8** in the final published article.

## Attribution note

This mapping identifies where code in this repository supports published panels. It does not imply sole authorship of those panels or of the complete computational workflow. The study was collaborative, and the official mouse scRNA-seq/TCR-seq pipeline is archived separately by Victor Gourain.
