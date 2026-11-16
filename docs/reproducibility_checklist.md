# Reproducibility checklist — UC-4

Reviewer request: *"To support reproducibility, the authors should deposit the
analysis code, software package versions, the sample inclusion list, annotations,
gene and probe mapping files, processed expression matrices, random seeds, and
the final model specification in a public repository."*

| # | Reviewer item | UC-4 location | Status |
|---|---|---|---|
| 1 | Analysis code | `00_`–`09_*.R`, `R/` (shared modules), `utils.R`, `run_all.R`; originals in `docs/legacy_scripts/` | ✓ |
| 2 | Software package versions | `config/package_versions.tsv` (116 packages), `sessionInfo.txt` (R 4.4.1 / Bioc 3.19) | ✓ |
| 3 | Sample inclusion list | `metadata/sample_inclusion_list.tsv` (written by `01_preprocess.R`) | ✓ |
| 4 | Annotations | `metadata/ANNOTATION_SOURCES.md`; annotation files under `data/raw/` | ✓ |
| 5 | Gene / probe mapping | `data/processed/<cohort>/expr_annotated_collapsed.tsv`; mapping logic in `R/01_preprocess.R` | ✓ |
| 6 | Processed expression matrices | `data/processed/<cohort>/expr_normalized.tsv` (one canonical copy per cohort) | ✓ |
| 7 | Random seeds | `config/seeds.yaml` (with origins); consumed via `SEED$...` | ✓ |
| 8 | Final model specification | `model/model_specification.md`, `model/final_model_spec.json`, `config/params.yaml` (`model:` section) | ✓ |

## Additional items requested by the author

| Item | UC-4 location |
|---|---|
| Deduplicated data (no repeated copies) | `metadata/file_manifest.tsv`, `docs/deduplication_report.md` |
| Merged code (same inputs → one script/function) | `R/` modules, `docs/analysis_workflow.md` §3 |
| Provenance of every UC-3 file | `metadata/file_manifest.tsv` |
| External web-tool parameters | `docs/external_web_tools.md` |
