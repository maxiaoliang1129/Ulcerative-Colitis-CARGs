# =============================================================================
#  04_risk_model.R
#  Fit the two-gene logistic model, full discrimination (ROC/AUC 95% CI,
#  single-gene AUC, Youden cutoff), nomogram, calibration and DCA.
#  Replaces blocks 1-5 of UC-3/06Model.R (01model, 02nom, 03dca) and the
#  separate per-cohort model.R copies.
#  Output: data/output/model/, model/final_model_spec.json, figures/model/
# =============================================================================
source(file.path(Sys.getenv("UC4_HOME", getwd()), "utils.R"))
suppressPackageStartupMessages(library(jsonlite))

expr <- read_expr(file.path(cohort_dir(TRAIN_COHORT), "expr_normalized.tsv"))
grp  <- read_group(file.path(cohort_dir(TRAIN_COHORT), "group.tsv"))

fit  <- fit_model(expr, grp)
prob <- apply_frozen_model(expr)
perf <- evaluate_model(prob, grp, TRAIN_COHORT)

write_tsv(data.frame(term = rownames(fit$coefficients), fit$coefficients, check.names = FALSE),
          file.path(DIR_TABLES, "model", "final_model_coefficients.tsv"))
write_tsv(data.frame(cohort = perf$cohort, n = perf$n, n_uc = perf$n_uc,
                     auc = perf$auc, ci_low = perf$ci[1], ci_high = perf$ci[3],
                     youden_cutoff = perf$cutoff),
          file.path(DIR_TABLES, "model", "model_performance_training.tsv"))

nomogram(fit, file.path(DIR_FIGURES, "model", "nomogram.pdf"))
cal <- calibration(fit)
dca(prob, as.integer(grp$Type == "UC"), file.path(DIR_FIGURES, "model", "dca_training.pdf"))

spec <- list(name = MODEL$name, formula = MODEL$formula, genes = MODEL$genes,
             coefficients = MODEL$coefficients, cutoff = MODEL$cutoff,
             r_version = R.version.string)
jsonlite::write_json(spec, file.path(DIR_MODEL, "final_model_spec.json"),
                     auto_unbox = TRUE, pretty = TRUE)
msg("training AUC = ", round(perf$auc, 4), " (", round(perf$ci[1], 4), "-", round(perf$ci[3], 4), ")")
