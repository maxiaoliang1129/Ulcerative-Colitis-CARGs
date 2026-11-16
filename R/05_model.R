# =============================================================================
#  R/05_model.R  --  risk model, ROC, nomogram, DCA (shared)
#
#  Consolidates the seven blocks of UC-3/06Model.R and the nine per-cohort
#  model.R copies.  One fit_model() and one evaluate_model() serve the
#  training cohort and all three validation cohorts.
# =============================================================================

model_matrix <- function(expr) {
  genes <- MODEL$genes
  x <- t(expr[intersect(genes, rownames(expr)), , drop = FALSE])
  if (!setequal(colnames(x), genes)) {
    missing <- setdiff(genes, rownames(expr))
    stop("Model genes missing from expression matrix: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  x[, genes, drop = FALSE]
}

fit_model <- function(expr, group) {
  x <- model_matrix(expr)
  grp <- factor(group$Type[match(rownames(x), group$Sample)], levels = c("control", "UC"))
  d <- data.frame(Group = as.integer(grp) - 1L, x, check.names = FALSE)
  suppressPackageStartupMessages(library(rms))
  fit <- glm(stats::as.formula(MODEL$formula), data = d, family = stats::binomial())
  cf <- stats::coef(summary(fit))
  list(fit = fit, coefficients = cf, data = d)
}

roc_auc <- function(prob, y, ci = TRUE) {
  suppressPackageStartupMessages(library(pROC))
  r <- pROC::roc(y, prob, quiet = TRUE)
  out <- list(auc = as.numeric(pROC::auc(r)), roc = r)
  if (ci) out$ci <- as.numeric(pROC::ci.auc(r))
  out
}

youden_cutoff <- function(prob, y) {
  r <- pROC::roc(y, prob, quiet = TRUE)
  coords <- pROC::coords(r, "best", best.method = "youden", ret = c("threshold","sensitivity","specificity"))
  as.numeric(coords["threshold"])
}

nomogram <- function(fit, path) {
  suppressPackageStartupMessages(library(rms))
  d <- fit$data
  dd <- datadist(d); options(datadist = "dd")
  lrm_fit <- rms::lrm(Group ~ CCL11 + CLDN1, data = d, x = TRUE, y = TRUE)
  nom <- rms::nomogram(lrm_fit, fun = function(x) 1 / (1 + exp(-x)))
  grDevices::pdf(path, width = 7, height = 7); on.exit(grDevices::dev.off()); plot(nom)
  invisible(lrm_fit)
}

calibration <- function(fit, B = MODEL$calibration_bootstraps) {
  suppressPackageStartupMessages(library(rms))
  d <- fit$data
  dd <- datadist(d); options(datadist = "dd")
  lrm_fit <- rms::lrm(Group ~ CCL11 + CLDN1, data = d, x = TRUE, y = TRUE)
  cal <- rms::calibrate(lrm_fit, method = "boot", B = B)
  list(fit = lrm_fit, calibration = cal)
}

dca <- function(prob, y, path = NULL) {
  suppressPackageStartupMessages(library(rmda))
  d <- data.frame(outcome = y, predictor = prob)
  dc <- rmda::decision_curve(formula = outcome ~ predictor, family = binomial(link = "logit"),
                             data = d, thresholds = seq(0, 1, by = 0.01), bootstraps = 50)
  if (!is.null(path)) {
    grDevices::pdf(path, width = 7, height = 6); on.exit(grDevices::dev.off())
    plot(dc)
  }
  dc
}

logistic_recalibration <- function(prob, y) {
  lp <- log(prob / (1 - prob))
  fit <- stats::glm(y ~ lp, family = stats::binomial())
  cf <- stats::coef(fit)
  recal <- plogis(cf[1] + cf[2] * lp)
  list(slope = unname(cf[2]), intercept = unname(cf[1]),
       prob = recal, brier = mean((recal - y)^2))
}

evaluate_model <- function(prob, group, cohort, out_dir = file.path(DIR_TABLES, "model")) {
  ensure_dir(out_dir)
  y <- as.integer(group$Type[match(names(prob), group$Sample)] == "UC")
  r <- roc_auc(prob, y)
  cut <- youden_cutoff(prob, y)
  list(cohort = cohort, auc = r$auc, ci = r$ci, cutoff = cut,
       roc = r$roc, n = length(y), n_uc = sum(y, na.rm = TRUE))
}

apply_frozen_model <- function(expr) {
  # Uses the published coefficients directly (model specification is frozen).
  x <- model_matrix(expr)
  lp <- as.numeric(MODEL$coefficients["(Intercept)"]) +
        as.matrix(x) %*% vapply(MODEL$genes, function(g) MODEL$coefficients[[g]], numeric(1))
  as.numeric(1 / (1 + exp(-lp)))
}
