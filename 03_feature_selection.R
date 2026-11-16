# =============================================================================
#  03_feature_selection.R
#  LASSO + random forest + XGBoost over one shared 1,320 x 35 design matrix,
#  then the three-way Venn/UpSet consensus.  Replaces 05ML.R and the separate
#  lasso.R / XGBoost.R / venn.R copies (identical inputs, redundant code).
#  Output: data/output/feature_selection/, figures/feature_selection/
# =============================================================================
source(file.path(Sys.getenv("UC4_HOME", getwd()), "utils.R"))

dat <- get_candidate_matrix()
msg("candidate matrix: ", nrow(dat$x), " samples x ", ncol(dat$x), " genes")

lasso <- run_lasso(dat)
rf    <- run_random_forest(dat)
xgb   <- run_xgboost(dat)

out <- file.path(DIR_TABLES, "feature_selection"); ensure_dir(out)
write_tsv(lasso$coef, file.path(out, "lasso_coefficients.tsv"))
writeLines(lasso$genes, file.path(out, "lasso_genes.txt"))
write_tsv(rf$importance, file.path(out, "rf_importance.tsv"))
writeLines(rf$genes, file.path(out, "rf_genes.txt"))
write_tsv(xgb$importance, file.path(out, "xgboost_importance.tsv"))
writeLines(xgb$genes, file.path(out, "xgboost_genes.txt"))

overlap <- Reduce(intersect, list(lasso$genes, rf$genes, xgb$genes))
writeLines(overlap, file.path(out, "ml_overlap_genes.txt"))
msg("three-way overlap: ", paste(overlap, collapse = ", "))
