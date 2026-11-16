# =============================================================================
#  R/06_immune.R  --  immune microenvironment (shared)
#
#  The ssGSEA -> ESTIMATE -> CIBERSORT -> correlation pipeline appears twice in
#  UC-3: once in 11Immue.R (UC vs control) and once in 12Subtype/03immue/
#  (subtype split, plus GEE).  This one parameterised module serves both:
#  immune_analysis(expr, grouping, ...) is called with group.txt in script 08
#  and with NMF_result.csv in script 07.
# =============================================================================

#' ssGSEA scores with GSVA against an immune gene set collection.
run_ssgsea <- function(expr, gmt = file.path(DIR_METADATA, "gene_sets", "immune.gmt")) {
  suppressPackageStartupMessages({ library(GSVA); library(GSEABase) })
  genesets <- GSEABase::getGmt(gmt)
  par <- GSVA::ssgseaParam(expr, genesets, normalize = TRUE)
  GSVA::gsvaParam # touch namespace
  GSVA::ssgsea(par)
}

#' ESTIMATE stromal / immune / estimate scores (vendored method).
run_estimate <- function(expr) {
  suppressPackageStartupMessages(library(estimate))
  common <- intersect(rownames(expr), estimate::commonGenes)
  x <- expr[common, , drop = FALSE]
  estimate::estimateScore(x, platform = "Affymetrix")
}

#' CIBERSORT deconvolution via the single vendored helper in lib/CIBERSORT.R.
run_cibersort <- function(expr, out_dir = file.path(DIR_TABLES, "immune")) {
  ensure_dir(out_dir)
  source(file.path(DIR_LIB, "CIBERSORT.R"), local = FALSE)
  sig <- read.delim(file.path(DIR_METADATA, "external", "cibersort_signature.txt"),
                    row.names = 1, check.names = FALSE)
  mix <- data.frame(Gene = rownames(expr), expr, check.names = FALSE)
  tmp_sig  <- file.path(out_dir, "signature.tmp.txt")
  tmp_mix  <- file.path(out_dir, "mixture.tmp.txt")
  write.table(sig, tmp_sig, sep = "\t", quote = FALSE)
  write.table(mix, tmp_mix, sep = "\t", quote = FALSE, row.names = FALSE)
  res <- CIBERSORT(tmp_sig, tmp_mix, perm = PARAMS$immune$cibersort$permutation)
  write_tsv(res, file.path(out_dir, "CIBERSORT-Results.tsv"))
  file.remove(tmp_sig, tmp_mix)
  res
}

#' Repeated-measures correlation between two genes and each immune fraction.
#' Requires a patient column; ignores patients with fewer than 2 observations.
run_rmcorr <- function(scores, patients, genes = MODEL$genes) {
  suppressPackageStartupMessages({ library(rmcorr); library(dplyr) })
  out <- list()
  min_obs <- PARAMS$immune$rmcorr$min_obs_per_patient
  for (g in genes) for (cell in setdiff(colnames(scores), g)) {
    d <- data.frame(Participant = patients, Measure1 = scores[[g]], Measure2 = scores[[cell]])
    d <- d[!is.na(d$Measure1) & !is.na(d$Measure2), ]
    n <- table(d$Participant)
    d <- d[d$Participant %in% names(n)[n >= min_obs], ]
    if (nrow(d) < 3) next
    fit <- tryCatch(rmcorr::rmcorr(participant = Participant, measure1 = Measure1,
                                   measure2 = Measure2, dataset = d), error = function(e) NULL)
    if (!is.null(fit)) out[[paste(g, cell, sep = "~")]] <-
      data.frame(gene = g, immune_cell = cell, r = fit$r, p = fit$p.value,
                 n_biopsies = nrow(d), n_patients = length(unique(d$Participant)))
  }
  res <- do.call(rbind, out)
  if (nrow(res)) res$q <- stats::p.adjust(res$p, "BH")
  res
}

#' One immune pipeline, parameterised by the grouping variable.
#' @param expr      genes x samples normalised matrix
#' @param grouping  factor named by sample (UC/control OR subtype)
#' @param patients  named patient ids (for rmcorr / GEE)
#' @param tag       output suffix, e.g. "group" or "subtype"
immune_analysis <- function(expr, grouping, patients = NULL,
                             tag = "group", out_dir = file.path(DIR_TABLES, "immune")) {
  ensure_dir(out_dir)
  msg("immune analysis [", tag, "]: ssGSEA")
  ss <- run_ssgsea(expr)
  write_tsv(data.frame(Sample = colnames(ss), t(ss), check.names = FALSE),
            file.path(out_dir, sprintf("ssgseaOut_%s.tsv", tag)))
  msg("immune analysis [", tag, "]: ESTIMATE")
  est <- run_estimate(expr)
  write_tsv(est, file.path(out_dir, sprintf("estimateScore_%s.tsv", tag)))

  msg("immune analysis [", tag, "]: CIBERSORT")
  cit <- tryCatch(run_cibersort(expr, out_dir), error = function(e) {
    msg("CIBERSORT failed: ", conditionMessage(e)); NULL })

  res <- list(ssgsea = ss, estimate = est, cibersort = cit)
  if (!is.null(cit) && !is.null(patients)) {
    common <- intersect(colnames(cit), colnames(ss))
    scores <- cbind(ss[, common, drop = FALSE], t(cit[, 1:min(22, ncol(cit)), drop = FALSE]))
    res$rmcorr <- run_rmcorr(as.data.frame(t(scores)), patients[colnames(ss)])
  }
  res
}

#' GEE-adjusted immune/ESTIMATE scores accounting for repeated biopsies.
#' Replaces three separate GEE blocks in 12Subtype.R.
run_gee <- function(scores, metadata, formula_rhs) {
  suppressPackageStartupMessages(library(geepack))
  d <- merge(metadata, scores, by = "Sample")
  out <- list()
  cols <- setdiff(colnames(scores), "Sample")
  for (col in cols) {
    f <- stats::as.formula(sprintf("%s ~ %s", col, formula_rhs))
    fit <- tryCatch(geepack::geeglm(f, data = d, id = d$Patient,
                                    family = stats::gaussian(),
                                    corstr = PARAMS$subtype$gee$corstr),
                    error = function(e) NULL)
    if (!is.null(fit)) out[[col]] <- data.frame(
      score = col, term = names(stats::coef(fit)), estimate = as.numeric(stats::coef(fit)),
      p.value = summary(fit)$coefficients[, "p.value"])
  }
  do.call(rbind, out)
}
