# =============================================================================
#  01_preprocess.R
#  Raw GEO files -> one canonical normalised matrix + group table per cohort.
#  All four cohorts are prepared through the single prepare_cohort() dispatcher
#  in R/01_preprocess.R (the four near-identical data*.R scripts of UC-3).
#  Output: data/processed/<cohort>/{expr_normalized,group,...}
# =============================================================================
source(file.path(Sys.getenv("UC4_HOME", getwd()), "utils.R"))
suppressPackageStartupMessages({ library(limma); library(edgeR); library(impute); library(data.table) })

incl <- list()
for (cohort in c(TRAIN_COHORT, VALIDATION_COHORTS)) {
  res <- prepare_cohort(cohort)
  if (!is.null(res$inclusion)) incl[[cohort]] <- res$inclusion
}

## one sample-inclusion list covering every cohort
if (length(incl)) {
  all_incl <- do.call(rbind, incl)
  write_tsv(all_incl[, intersect(c("cohort","Sample","Patient","Type","included","exclusion_reason"),
                                  colnames(all_incl))],
            file.path(DIR_METADATA, "sample_inclusion_list.tsv"))
}
msg("preprocessing complete: ", length(COHORTS), " cohorts")
