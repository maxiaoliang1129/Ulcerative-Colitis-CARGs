# =============================================================================
#  R/01_preprocess.R  --  cohort preparation functions (shared)
#
#  Consolidates the four near-identical preparation scripts in UC-3 into three
#  parameterised routes:
#    qc_rnaseq_ncbi()      <- 00data/00train/GSE193677/data.R
#    qc_rnaseq_deseq2()    <- 00data/01test/GSE128682/00data/data1.R
#                             (byte-identical copy in 06Model/04verify/...)
#    qc_affymetrix()       <- GSE87466/data1.R and GSE206285/data.R (same logic)
#  QC boxplots, raw/collapsed matrices, group tables and the one-biopsy-per-
#  patient sampling are produced by the same functions for every cohort.
#  Source R/00_config.R first.
# =============================================================================

run_qc_plots <- function(expr, cohort, prefix) {
  ensure_dir(file.path(UC4_ROOT, "figures", "qc"))
  df <- as.data.frame(expr)
  p <- ggplot2::ggplot(reshape2::melt(df), ggplot2::aes_string(x = "variable", y = "value")) +
    ggplot2::geom_boxplot(ggplot2::aes_string(group = "variable"), fill = "steelblue", alpha = 0.7) +
    ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Sample", y = "Expression", title = sprintf("%s %s", cohort, prefix))
  save_pdf(p, file.path(UC4_ROOT, "figures", "qc", sprintf("%s_%s.pdf", prefix, cohort)), 12, 5)
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## RNA-seq counts with NCBI GRCh38 annotation (training cohort)
## ---------------------------------------------------------------------------
qc_rnaseq_ncbi <- function(cohort = TRAIN_COHORT) {
  cfg <- COHORTS[[cohort]]
  raw <- file.path(DIR_RAW, cohort)
  out <- cohort_dir(cohort)
  msg("preparing ", cohort, " (rnaseq_ncbi)")

  counts <- data.table::fread(require_raw(cohort, cfg$files$counts), sep = "\t",
                              header = TRUE, data.table = FALSE, check.names = FALSE)
  colnames(counts)[1] <- "ID"
  annot <- read_annotation(require_raw(cohort, cfg$files$annotation),
                           id_col = "GeneID", symbol_col = "Symbol")
  expr <- as.matrix(counts[, -1, drop = FALSE]); rownames(expr) <- counts$ID
  storage.mode(expr) <- "double"
  if (isTRUE(cfg$qc_log2)) expr <- log2(expr + PREP$log2_pseudocount)
  run_qc_plots(expr, cohort, "rawBox")

  expr <- map_probes_to_symbols(expr, annot)
  expr_collapsed <- collapse_duplicate_symbols(expr)
  write_tsv(data.frame(Gene = rownames(expr_collapsed), expr_collapsed, check.names = FALSE),
            file.path(out, "expr_annotated_collapsed.tsv"), row.names = FALSE)

  ## phenotype and patient bookkeeping
  suppressPackageStartupMessages(library(GEOquery))
  gse <- GEOquery::getGEO(filename = require_raw(cohort, cfg$files$phenotype), getGPL = FALSE)
  pheno <- Biobase::phenoData(gse)@data
  pheno$patient_id <- sub(".*reGRID_([0-9]+)_Biopsy_.*", "\\1", pheno$title)

  col <- cfg$group_column
  lab <- pheno[[col]]; lab[lab == "Control"] <- "control"
  grp <- data.frame(Sample = pheno$geo_accession, Type = lab, stringsAsFactors = FALSE)
  grp <- grp[grp$Type %in% c("UC", "control"), ]
  aln <- align_group(expr_collapsed, grp)
  write_tsv(aln$group, file.path(out, "group.tsv"))

  ## one biopsy per patient (seeded)
  sp <- data.frame(Sample = pheno$geo_accession, Patient = pheno$patient_id,
                   Type = lab, stringsAsFactors = FALSE)
  write_tsv(sp, file.path(out, "sample_patient_group.tsv"))
  set.seed(SEED$one_biopsy_per_patient)
  one <- do.call(rbind, lapply(split(sp, sp$Patient), function(d) d[sample(nrow(d), 1), ]))
  write_tsv(one, file.path(out, "group_one_biopsy_per_patient.tsv"))
  write_tsv(grp[grp$Type == "UC", ], file.path(out, "group_all_samples.tsv"))

  ## sample inclusion list row for this cohort
  incl <- transform(sp, cohort = cohort, exclusion_reason = NA_character_)
  incl$included <- incl$Sample %in% aln$group$Sample

  norm <- normalise_expression(aln$expr, log2 = TRUE, impute = TRUE)
  write_tsv(data.frame(Gene = rownames(norm), norm, check.names = FALSE),
            file.path(out, "expr_normalized.tsv"))
  run_qc_plots(norm, cohort, "normalBox")
  write_tsv(qc_stats(norm), file.path(DIR_TABLES, "qc", sprintf("qc_%s.tsv", cohort)))
  list(expr = norm, group = aln$group, inclusion = incl)
}

## ---------------------------------------------------------------------------
## DESeq2-normalised counts (GSE128682).  Labels come from the sample-name
## prefix, which fixes the positional mis-assignment in the original script.
## ---------------------------------------------------------------------------
qc_rnaseq_deseq2 <- function(cohort = "GSE128682") {
  cfg <- COHORTS[[cohort]]
  out <- cohort_dir(cohort)
  msg("preparing ", cohort, " (rnaseq_deseq2)")
  counts <- data.table::fread(require_raw(cohort, cfg$files$counts), sep = "\t",
                              header = TRUE, data.table = FALSE, check.names = FALSE)
  colnames(counts)[1] <- "ID"
  expr <- as.matrix(counts[, -1, drop = FALSE]); rownames(expr) <- counts$ID
  storage.mode(expr) <- "double"
  expr <- log2(expr + PREP$log2_pseudocount)
  run_qc_plots(expr, cohort, "rawBox")

  expr <- expr[!duplicated(toupper(rownames(expr))), , drop = FALSE]
  rownames(expr) <- toupper(rownames(expr))
  write_tsv(data.frame(Gene = rownames(expr), expr, check.names = FALSE),
            file.path(out, "expr_annotated_collapsed.tsv"))

  lab <- ifelse(grepl("^N_", colnames(expr)), "control",
                ifelse(grepl("^UC_", colnames(expr)), "UC", NA_character_))
  grp <- data.frame(Sample = colnames(expr), Type = lab, stringsAsFactors = FALSE)
  grp <- grp[!is.na(grp$Type), ]
  aln <- align_group(expr, grp)
  write_tsv(aln$group, file.path(out, "group.tsv"))

  norm <- normalise_expression(aln$expr, log2 = FALSE, impute = FALSE)
  write_tsv(data.frame(Gene = rownames(norm), norm, check.names = FALSE),
            file.path(out, "expr_normalized.tsv"))
  run_qc_plots(norm, cohort, "normalBox")
  list(expr = norm, group = aln$group)
}

## ---------------------------------------------------------------------------
## Affymetrix series matrices (GSE87466 and GSE206285) -- one shared route
## ---------------------------------------------------------------------------
qc_affymetrix <- function(cohort) {
  cfg <- COHORTS[[cohort]]
  out <- cohort_dir(cohort)
  msg("preparing ", cohort, " (affymetrix)")
  eset <- GEOquery::getGEO(filename = require_raw(cohort, cfg$files$expression),
                           getGPL = FALSE, parsingfeat = FALSE)
  expr <- Biobase::exprs(eset)
  annot <- read_annotation(file.path(DIR_RAW, "GPL13158", "GPL13158-5065.txt"),
                           id_col = "ID", symbol_col = "Gene symbol")
  expr <- map_probes_to_symbols(expr, annot)
  expr_collapsed <- collapse_duplicate_symbols(expr)
  write_tsv(data.frame(Gene = rownames(expr_collapsed), expr_collapsed, check.names = FALSE),
            file.path(out, "expr_annotated_collapsed.tsv"))

  pheno <- Biobase::phenoData(eset)@data
  col <- cfg$group_column
  lab <- pheno[[col]]
  if (!is.null(cfg$group_map)) for (k in names(cfg$group_map)) lab[lab == k] <- cfg$group_map[[k]]
  grp <- data.frame(Sample = pheno$geo_accession, Type = lab, stringsAsFactors = FALSE)
  grp <- grp[grp$Type %in% c("UC", "control"), ]
  aln <- align_group(expr_collapsed, grp)
  write_tsv(aln$group, file.path(out, "group.tsv"))

  # impute.knn is applied identically to both cohorts (fixes C2 in UC-2)
  norm <- normalise_expression(aln$expr, log2 = FALSE, impute = TRUE)
  write_tsv(data.frame(Gene = rownames(norm), norm, check.names = FALSE),
            file.path(out, "expr_normalized.tsv"))
  run_qc_plots(norm, cohort, "normalBox")
  list(expr = norm, group = aln$group)
}

PREP_ROUTES <- list(rnaseq_ncbi = qc_rnaseq_ncbi,
                    rnaseq_deseq2 = qc_rnaseq_deseq2,
                    affymetrix = qc_affymetrix)

prepare_cohort <- function(cohort) {
  route <- COHORTS[[cohort]]$route
  if (is.null(route)) stop("No route configured for ", cohort, call. = FALSE)
  PREP_ROUTES[[route]](cohort)
}
