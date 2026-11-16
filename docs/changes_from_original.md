# Changes relative to the original analysis (UC-3) — UC-4

This document records what was removed, merged or corrected when UC-3 was
restructured into UC-4. See also `deduplication_report.md` (file-level MD5
dedup) and `analysis_workflow.md` §3 (code consolidation map).

## 1. Removed because byte-identical (MD5)

361 source files → 228 unique files; 133 copies removed (13.56 GB).

| File | Copies in UC-3 | Kept in UC-4 |
|---|---:|---|
| training `normalize.txt` (924 MB) | 12 | 1 → `data/processed/GSE193677/expr_normalized.txt` |
| training `group.txt` | 17 | 1 → `data/processed/GSE193677/group.tsv` |
| 2-gene `final.model.txt` / `venn.txt` | 16 | 1 → `data/derived/model/final_model_genes.txt` |
| `diffGeneExp_p0.05.txt` (DEG matrix) | 9 | 1 → `data/derived/deg/diffGeneExp_p0.05.txt` |
| `NMF_result.csv` (subtype calls) | 7 | 1 → `data/derived/subtype/NMF_result.csv` |
| `GPL13158-5065.txt` (58 MB) | 3 | 1 → `data/raw/GPL13158/` |
| `ssGSEA18.CIBERSORT.R` | 2 | 1 → `lib/CIBERSORT.R` |
| per-cohort `normalize/group` in `06Model/04verify/` | 3–5 each | 1 each |

## 2. Merged because inputs were identical

| Original | Consolidated into |
|---|---|
| 4 × `data.R` / `data1.R` | `R/01_preprocess.R` (three routes + one dispatcher) |
| 3 DEG blocks in `01.limma_diff.R` | `R/02_diff.R::run_limma()` (+ block argument) |
| 4 GO/KEGG blocks (`03.venn.R`, `04String.R`, `12Subtype.R`) | `R/03_enrichment.R::enrichment()` |
| `05ML.R` + `lasso.R` + `XGBoost.R` | `R/04_featuresel.R` (one shared design matrix) |
| 7 blocks of `06Model.R` + 9 × `model.R` / `boxplot.R` | `R/05_model.R`, scripts 04–05 (cohort loop) |
| `11Immue.R` + `estimate.R` + `12Subtype/03immue/` immune code | `R/06_immune.R` (one pipeline, grouping argument) |

The ten original stage scripts and all nested helpers are kept verbatim in
`docs/legacy_scripts/`.

## 3. Corrections applied in the consolidated code

| # | Problem in UC-3 | UC-4 behaviour |
|---|---|---|
| C1 | GSE128682 labels assigned by position against differently ordered columns | derived from sample-name prefix (`^N_` → control, `^UC_` → UC) |
| C2 | GSE87466 `impute.knn()` result immediately overwritten | both Affymetrix cohorts imputed identically |
| C3 | `path <- 'D:\\LZ-25234\\UC\\...'; setwd()` in every script | paths resolved from project root / `config/params.yaml` |
| C4 | 89 invalid RStudio `'''` block delimiters | removed; all R files parse under R 4.4.1 |
| C5 | wrong-project path `D:\LZ-25197\sepsis\...` in `04String.R` | removed |
| C6 | two different gene lists both named `diff_intergene.txt` | disambiguated by canonical destination |
| C7 | `height = 3。5` (full-width period) in `06Model.R` | removed with the archived script |
| C8 | package versions not recorded | `config/package_versions.tsv` + `sessionInfo.txt` |

## 4. Reproducibility notes

* All randomness from `config/seeds.yaml`; no literal seed in active code.
* Every original file is traceable through `metadata/file_manifest.tsv`.
* Two version-sensitive numbers (duplicateCorrelation consensus; 8/9 feature
  genes) are documented in `analysis_workflow.md` §6; the two model genes and
  all reported AUCs are unaffected.
