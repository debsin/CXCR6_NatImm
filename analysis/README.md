# Curated analysis code

This directory contains the publication-relevant computational analyses from this repository, reorganized by analytical objective. The code itself is preserved from the original 2024 repository; curation here changes presentation and scope, not the underlying scientific logic.

## Analysis blocks

| Directory | Published role | Main data source(s) | Reproducibility status |
| --- | --- | --- | --- |
| `01_bal_covid_validation/` | Figure 8a–b; Extended Data Fig. 10a | GSE145926 | Public-data-derived analysis, but the preserved script expects locally prepared Seurat objects |
| `02_bcg_trained_immunity/` | Figure 8c | GSE124218 / GSE124220 | Public-data-derived analysis; reconstruction of the mouse AM signature is required |
| `03_cross_species_trm/` | Figure 8d–e; Extended Data Fig. 10e–f | Mouse scRNA-seq, public TRM signatures, EGAS00001004707 | Requires study-specific mouse objects and, for human tumor TRM data, controlled-access acquisition |
| `04_luad_chemokine_correlations/` | Figure 8f; Extended Data Fig. 10d | TCGA LUAD / UCSC Xena; published TRM signatures | Public-data-derived; paths and inputs must be configured locally |
| `05_pan_cancer_survival/` | Figure 8g–j | TCGA / UCSC Xena | Public-data-derived; paths and inputs must be configured locally |

## Preserved filenames

Several script names retain manuscript-development numbering such as `NatImm_Figure9.R`. The corresponding human-validation figure became **Figure 8** in the final Nature Immunology article.

## Important reproducibility note

These scripts originated within an active collaborative analysis environment. Some contain absolute `C:/Projects/...` paths, assume objects created by earlier preprocessing, or use package interfaces current at the time of analysis. They are therefore presented as transparent scientific analysis code, not as a one-command software package.

The repository intentionally does **not** include patient-level data, processed RDS objects, TCGA matrices, controlled-access files or publisher source-data spreadsheets. See [`../data/README.md`](../data/README.md) and [`../docs/data-sources.md`](../docs/data-sources.md).

The complete original working layout is preserved on `archive/original-2024`.
