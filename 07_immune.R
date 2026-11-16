# =============================================================================
#  07_immune.R
#  Immune microenvironment of the UC-vs-control comparison: ssGSEA, ESTIMATE,
#  CIBERSORT, immune-boxplots and rmcorr / patient-level Spearman correlations.
#  Uses the SAME immune_analysis() function as 06_subtype.R (the second copy of
#  this pipeline in UC-3/12Subtype/03immue) -- only the grouping differs.
#  Replaces 11Immue.R (738 lines) and the estimate.R duplicate.
#  Output: data/output/immune/, figures/immune/
# =============================================================================
source(file.path(Sys.getenv("UC4_HOME", getwd()), "utils.R"))
out <- ensure_dir(file.path(DIR_TABLES, "immune"))

expr <- read_expr(file.path(cohort_dir(TRAIN_COHORT), "expr_normalized.tsv"))
grp  <- read_group(file.path(cohort_dir(TRAIN_COHORT), "group.tsv"))
sp   <- read_group(file.path(cohort_dir(TRAIN_COHORT), "sample_patient_group.tsv"))
patients <- stats::setNames(sp$Patient, sp$Sample)

imm <- immune_analysis(expr, grouping = stats::setNames(grp$Type, grp$Sample),
                       patients = patients[colnames(expr)], tag = "group", out_dir = out)
if (!is.null(imm$rmcorr)) {
  write_tsv(imm$rmcorr, file.path(out, "rmcorr_gene_immune.tsv"))
  ## patient-level Spearman as a complementary summary
  pl <- stats::aggregate(as.data.frame(t(imm$cibersort[, 1:22, drop = FALSE])),
                         list(Patient = patients[colnames(imm$cibersort)]), mean)
  write_tsv(pl, file.path(out, "patient_level_immune_means.tsv"))
}
msg("immune analysis complete")
