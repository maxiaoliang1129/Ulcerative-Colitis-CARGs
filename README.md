# UC-4 — Transcriptomic Analysis of Ulcerative Colitis

## 1. Project Description

This project investigates transcriptomic diagnostic biomarkers for ulcerative
colitis (UC) and identifies a two-gene classifier, **CCL11 + CLDN1**:

```
logit P(UC) = -3.900392 + 0.417981 × CCL11 + 0.409387 × CLDN1
```

The model was trained in GSE193677 (1,320 biopsies from 660 patients; decision
threshold 0.639643) and validated in three external cohorts: GSE128682,
GSE87466 and GSE206285. The analyses cover differential expression, feature
selection, risk modelling and validation, molecular subtyping, immune
microenvironment characterization, and downstream ceRNA / TF / drug network
analyses. See `model/model_specification.md` for the full model specification.

## 2. Directory Structure

```
UC-4/
├── README.md                    this file
├── sessionInfo.txt              R session information with package versions
│
├── 00_download_data.R           download GEO datasets and platform annotations
├── 01_preprocess.R              raw data to normalized matrices and group tables
├── 02_diff_expression.R         differential expression, sensitivity and enrichment
├── 03_feature_selection.R       LASSO / random forest / XGBoost feature selection
├── 04_risk_model.R              two-gene risk model, ROC, nomogram, calibration, DCA
├── 05_validation.R              validation in external cohorts
├── 06_subtype.R                 consensus clustering and subtype analyses
├── 07_immune.R                  ssGSEA / ESTIMATE / CIBERSORT / correlation analyses
├── 08_downstream.R              single-gene GSEA, TF / ceRNA / drug networks
├── 09_figures.R                 combined figures and SourceData workbooks
│
├── R/                           shared function modules
│   ├── 00_config.R              path resolution, configuration and I/O helpers
│   ├── 01_preprocess.R          cohort preprocessing functions
│   ├── 02_diff.R                differential expression functions
│   ├── 03_enrichment.R          GO / KEGG enrichment functions
│   ├── 04_featuresel.R          feature selection functions
│   ├── 05_model.R               model fitting and evaluation functions
│   └── 06_immune.R              immune analysis functions
│
├── data/                        data files (157)
│   ├── raw/                     raw GEO datasets and platform files (14)
│   ├── processed/               normalized matrices and group tables (21)
│   └── derived/                 analysis result tables (122), grouped by analysis:
│                                deg / enrichment / feature_selection / model /
│                                immune / subtype / ppi / cerna / tf / drugs /
│                                gene_sets / organized / supplementary / other / logs
│
├── metadata/                    metadata and annotations (10 files)
│   ├── file_manifest.tsv        data source and file path record
│   ├── ANNOTATION_SOURCES.md    sources of annotations and gene/probe mappings
│   ├── gene_sets/               gene sets (ARGs, DEGs, immune.gmt, etc.)
│   └── external/                external reference data (CIBERSORT signatures,
│                                miRTarBase, etc.)
│
├── model/
│   └── model_specification.md   final model specification
│
├── lib/
│   ├── CIBERSORT.R              CIBERSORT implementation
│   └── ConsensusClusterPlus_fixed.R   consensus clustering implementation
│
└── docs/                        documentation (30 files)
    ├── analysis_workflow.md         analysis workflow
    ├── reproducibility_checklist.md reproducibility checklist
    ├── deduplication_report.md      data processing report
    ├── changes_from_original.md     change log
    ├── external_web_tools.md        external web tools
    └── legacy_scripts/              archived original R scripts
```

## 3. Software Environment

- **R version**: R 4.4.1 (2024-06-14 ucrt)
- **Platform**: x86_64-w64-mingw32/x64, Windows 11 x64
- **Bioconductor**: 3.19
- **R packages**: the complete list of packages and their exact versions used in
  this analysis is recorded in `sessionInfo.txt`.
