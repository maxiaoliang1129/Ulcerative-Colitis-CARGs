# =============================================================================
#  R/02_diff.R  --  differential expression (shared)
#
#  One run_limma() function replaces the three near-identical DEG blocks in
#  UC-3/01.limma_diff.R (main analysis, duplicateCorrelation sensitivity and
#  one-biopsy-per-patient sensitivity).  The block argument switches the
#  patient-correlation design without duplicating any fitting code.
# =============================================================================

#' Fit limma on an expression matrix with a binary group factor.
#' @param expr          genes x samples numeric matrix (normalised, log2)
#' @param group         data.frame with columns Sample, Type ("UC"/"control")
#' @param block         optional patient factor (for duplicateCorrelation)
run_limma <- function(expr, group, block = NULL) {
  common <- intersect(colnames(expr), group$Sample)
  expr <- expr[, common, drop = FALSE]
  grp <- factor(group$Type[match(common, group$Sample)], levels = c("UC", "control"))
  design <- stats::model.matrix(~ grp)
  colnames(design) <- c("intercept", "UC_vs_control")

  fit <- if (!is.null(block)) {
    b <- block[match(common, names(block))]
    corfit <- limma::duplicateCorrelation(expr, design, block = b)
    limma::lmFit(expr, design, block = b, correlation = corfit$consensus)
  } else {
    limma::lmFit(expr, design)
  }
  fit <- limma::eBayes(fit)
  list(fit = fit, design = design)
}

#' Assemble the full DEG table, then apply the configured thresholds.
degs_from_fit <- function(fit, coef = "UC_vs_control") {
  tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "P", adjust.method = "BH")
  tt <- tt[order(tt$adj.P.Val, -abs(tt$logFC)), , drop = FALSE]
  sig <- tt[tt$adj.P.Val < THRESH$deg_adj_p & abs(tt$logFC) > THRESH$deg_logfc, , drop = FALSE]
  list(all = tt, significant = sig)
}

write_deg_results <- function(res, prefix, out_dir = file.path(DIR_TABLES, "deg")) {
  ensure_dir(out_dir)
  write_tsv(data.frame(Gene = rownames(res$all), res$all, check.names = FALSE),
            file.path(out_dir, sprintf("%s_all_genes.tsv", prefix)), row.names = FALSE)
  write_tsv(data.frame(Gene = rownames(res$significant), res$significant, check.names = FALSE),
            file.path(out_dir, sprintf("%s.tsv", prefix)), row.names = FALSE)
  invisible(res)
}

volcano <- function(res, title, path) {
  d <- data.frame(logFC = res$all$logFC,
                  negLogP = -log10(res$all$adj.P.Val + 1e-300),
                  sig = rownames(res$all) %in% rownames(res$significant))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = logFC, y = negLogP, colour = sig)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.7) +
    ggplot2::scale_colour_manual(values = c("grey70", "red"), guide = FALSE) +
    ggplot2::theme_bw() + ggplot2::labs(title = title, x = "log2 fold change",
                                        y = "-log10(adjusted P)")
  save_pdf(p, path, 7, 6)
}

## ---- one shared function for the three designs ----------------------------
run_deg_analysis <- function(expr, group, patients = NULL, out_dir = file.path(DIR_TABLES, "deg")) {
  ensure_dir(out_dir)
  ## (A) main analysis
  main <- degs_from_fit(run_limma(expr, group)$fit)
  write_deg_results(main, "DEG_main", out_dir)
  volcano(main, "Main UC vs control", file.path(DIR_FIGURES, "deg", "Volcano_main.pdf"))

  ## (B) duplicateCorrelation sensitivity (blocked on patient)
  results <- list(main = main)
  if (!is.null(patients)) {
    b <- stats::setNames(factor(patients), colnames(expr))
    dup <- degs_from_fit(run_limma(expr, group, block = b)$fit)
    write_deg_results(dup, "DEG_dupcor", out_dir)
    volcano(dup, "duplicateCorrelation", file.path(DIR_FIGURES, "deg", "Volcano_dupcor.pdf"))
    results$dupcor <- dup
  }

  ## (C) one biopsy per patient (seeded subset)
  sp <- read_group(file.path(cohort_dir(TRAIN_COHORT), "sample_patient_group.tsv"))
  if (nrow(sp)) {
    set.seed(SEED$one_biopsy_per_patient)
    keep <- unname(tapply(sp$Sample, sp$Patient, function(z) sample(z, 1)))
    sub <- run_limma(expr[, intersect(colnames(expr), keep), drop = FALSE], group)
    ob <- degs_from_fit(sub$fit)
    write_deg_results(ob, "DEG_one_biopsy_per_patient", out_dir)
    volcano(ob, "One biopsy per patient", file.path(DIR_FIGURES, "deg", "Volcano_onePerPatient.pdf"))
    results$one_biopsy <- ob
  }
  results
}
