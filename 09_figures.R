# =============================================================================
#  09_figures.R
#  Regenerate every figure from the saved result tables (never from scratch),
#  and write one SourceData workbook per figure.  Replaces the plotting blocks
#  scattered through all ten UC-3 scripts; all plotting helpers live in R/.
#  Figures produced by earlier steps (Volcano, enrichment dotplots, CIBERSORT,
#  consensus clustering) are already in figures/<topic>/ -- this script adds
#  the combined model panels and assembles the source-data workbooks.
#  Output: figures/fig*.pdf, figures/source_data/
# =============================================================================
source(file.path(Sys.getenv("UC4_HOME", getwd()), "utils.R"))
suppressPackageStartupMessages(library(openxlsx))
ensure_dir(file.path(DIR_FIGURES, "source_data"))

## ---- combined ROC panel across all cohorts --------------------------------
roc_frames <- list()
paths <- c(file.path(DIR_TABLES, "model", "model_performance_training.tsv"),
           file.path(DIR_TABLES, "validation", "model_AUC_all_cohorts.tsv"))
for (p in paths) if (file.exists(p)) roc_frames[[length(roc_frames)+1]] <- read.delim(p)
perf <- if (length(roc_frames)) do.call(rbind, roc_frames) else NULL

if (!is.null(perf) && "auc" %in% colnames(perf)) {
  p <- ggplot2::ggplot(perf, ggplot2::aes(x = reorder(cohort, auc), y = auc)) +
    ggplot2::geom_col(fill = "steelblue") + ggplot2::coord_flip() +
    ggplot2::ylim(0.5, 1) + ggplot2::theme_bw() +
    ggplot2::labs(title = "Discrimination by cohort", x = NULL, y = "AUC")
  save_pdf(p, file.path(DIR_FIGURES, "fig_ROC_summary.pdf"))
}

## ---- source-data workbook per result table --------------------------------
tables <- list.files(file.path(DIR_TABLES), recursive = TRUE, pattern = "\\.tsv$", full.names = TRUE)
for (i in seq_along(tables)) {
  d <- tryCatch(read.delim(tables[i]), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) next
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "data")
  openxlsx::writeData(wb, "data", d)
  openxlsx::saveWorkbook(wb, file.path(DIR_FIGURES, "source_data",
                         sprintf("SourceData_%02d.xlsx", i)), overwrite = TRUE)
}
msg("figures and source-data workbooks written to ", DIR_FIGURES)
