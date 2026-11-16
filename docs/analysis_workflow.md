# Analysis workflow — UC-4

## 1. Requirements

* R 4.4.1 and Bioconductor 3.19 (recorded in `config/package_versions.tsv` and `sessionInfo.txt`).
* 16 GB RAM recommended (the training matrix is 39,374 × 1,320).
* 4+ CPU cores; only steps 02/06 use parallel workers.
* ~8 GB disk for the de-duplicated data plus generated outputs.

Raw GEO files are expected under `data/raw/`; the shared GPL13158 platform table is
stored **once** at `data/raw/GPL13158/GPL13158-5065.txt`. Set `UC4_RAW` to use an
existing raw-data directory instead of downloading.

## 2. Run order

| Step | Script | Reads | Writes | Notes |
|---|---|---|---|---|
| 00 | `00_download_data.R` | GEO (internet) | `data/raw/` | skip if raw files present |
| 01 | `01_preprocess.R` | `data/raw/` | `data/processed/<cohort>/` | all four cohorts via one dispatcher |
| 02 | `02_diff_expression.R` | training matrix | `data/output/deg/`, `gene_sets/`, `enrichment/` | main + 2 sensitivity DEGs |
| 03 | `03_feature_selection.R` | DEG + gene sets | `data/output/feature_selection/` | LASSO / RF / XGBoost |
| 04 | `04_risk_model.R` | training matrix | `data/output/model/`, `model/final_model_spec.json` | ROC, nomogram, DCA |
| 05 | `05_validation.R` | validation matrices | `data/output/validation/` | loop over 3 cohorts + recalibration |
| 06 | `06_subtype.R` | training matrix | `data/output/subtype/` | consensus clustering + subtype immune/GEE |
| 07 | `07_immune.R` | training matrix | `data/output/immune/` | ssGSEA, ESTIMATE, CIBERSORT, rmcorr |
| 08 | `08_downstream.R` | training matrix | `data/output/downstream/` | single-gene GSEA + web-tool registration |
| 09 | `09_figures.R` | all tables | `figures/` | combined figures + SourceData workbooks |

Run everything with `Rscript run_all.R`; run selected steps with
`Rscript run_all.R 02 03 04`. Each step is a separate `Rscript --vanilla`
process; logs are in `data/output/logs/`.

## 3. Code consolidation map

The original 28 R scripts collapse into 10 numbered scripts + 6 shared modules:

| Original UC-3 scripts | UC-4 location |
|---|---|
| `data.R`, `data1.R` (×4, near-identical) | `R/01_preprocess.R` (3 routes) |
| `01.limma_diff.R` (3 DEG blocks) | `R/02_diff.R::run_limma()` |
| `02.genesets.R`, `03.venn.R`, `04String.R` (4 GO/KEGG blocks) | `R/03_enrichment.R::enrichment()` |
| `05ML.R` + `lasso.R` + `XGBoost.R` | `R/04_featuresel.R` |
| `06Model.R` (7 blocks) + 9 × `model.R` | `R/05_model.R`, scripts 04–05 |
| `10ssgsea.R` | `08_downstream.R` |
| `11Immue.R` + `estimate.R` + `12Subtype/03immue/` | `R/06_immune.R`, scripts 06–07 |
| `12Subtype.R` (2,204 lines) | `06_subtype.R` + shared modules |

Original scripts are archived verbatim under `docs/legacy_scripts/`.

## 4. Randomness

Every stochastic step takes its seed from `config/seeds.yaml` via `SEED$...`;
no script contains a literal seed. Origins (original UC-3 line numbers) are
recorded as comments in that file.

## 5. What is not reproduced here

The TF network (`07TFs`), ceRNA network (`08ceRNA`) and drug repurposing
(`09DsigDB`) were produced with external web tools (NetworkAnalyst, miRTarBase,
DSigDB) and are shipped as their result files under `data/derived/{tf,cerna,drugs}/`,
with the tool parameters documented in `docs/external_web_tools.md`. WGCNA hub
genes were a fixed input (`metadata/gene_sets/wgcna_hubs.txt`); no WGCNA script
existed in the original tree.

## 6. Known version-sensitive results

Earlier verification on R 4.4.1 showed two numbers shift slightly with package
version: the `duplicateCorrelation` consensus (0.1503 vs 0.1422) and the
three-method feature intersection (8 of the 9 published genes; MMP12 vs CCN4).
The two model genes (CCL11, CLDN1) and all reported AUCs are unaffected.
