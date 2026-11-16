# =============================================================================
#  05_validation.R
#  Apply the frozen model to the three external cohorts in one loop, replacing
#  the nine near-identical model.R files under 00data/01test/* and
#  06Model/04verify/*.  Includes logistic recalibration (UC-3 06model2).
#  Output: data/output/validation/
# =============================================================================
source(file.path(Sys.getenv("UC4_HOME", getwd()), "utils.R"))
out <- ensure_dir(file.path(DIR_TABLES, "validation"))
rows <- list()
for (cohort in VALIDATION_COHORTS) {
  expr <- read_expr(file.path(cohort_dir(cohort), "expr_normalized.tsv"))
  grp  <- read_group(file.path(cohort_dir(cohort), "group.tsv"))
  prob <- apply_frozen_model(expr)
  perf <- evaluate_model(prob, grp, cohort)
  y <- as.integer(grp$Type[match(names(prob), grp$Sample)] == "UC")
  rc <- logistic_recalibration(prob, y)
  rows[[cohort]] <- data.frame(cohort = cohort, n = perf$n, n_uc = perf$n_uc,
    auc = perf$auc, ci_low = perf$ci[1], ci_high = perf$ci[3],
    brier = mean((prob - y)^2), recal_slope = rc$slope,
    recal_intercept = rc$intercept, brier_recalibrated = rc$brier,
    stringsAsFactors = FALSE)
  write_tsv(data.frame(sample = names(prob), group = grp$Type[match(names(prob), grp$Sample)],
                       probability = prob, check.names = FALSE),
            file.path(out, sprintf("predictions_%s.tsv", cohort)))
}
write_tsv(do.call(rbind, rows), file.path(out, "model_AUC_all_cohorts.tsv"))
msg("validation complete: ", paste(VALIDATION_COHORTS, collapse = ", "))
