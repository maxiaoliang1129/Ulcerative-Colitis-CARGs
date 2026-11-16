# =============================================================================
#  R/00_config.R  --  shared configuration, paths and I/O helpers
#
#  This single file replaces the hard-coded path blocks that appeared at the
#  top of every one of the ten UC-3 analysis scripts (D:\LZ-25234\UC\...,
#  E:\LZ-25007\colitis\..., setwd(), and the wrong-project sepsis path in
#  04String.R).  Every script sources this file and resolves all paths from
#  the project root.  Set UC4_HOME to override the auto-detected root.
# =============================================================================

.find_root <- function(p = getwd()) {
  for (i in seq_len(10)) {
    if (file.exists(file.path(p, "config", "params.yaml"))) return(p)
    if (identical(dirname(p), p)) break
    p <- dirname(p)
  }
  stop("Could not find the UC-4 project root. Run from inside the repository ",
       "or set UC4_HOME.", call. = FALSE)
}

UC4_ROOT <- Sys.getenv("UC4_HOME", unset = .find_root())
CONF_DIR <- file.path(UC4_ROOT, "config")

if (!requireNamespace("yaml", quietly = TRUE))
  stop("Package 'yaml' is required to read the configuration.", call. = FALSE)

PARAMS <- yaml::yaml.load_file(file.path(CONF_DIR, "params.yaml"))
SEED   <- yaml::yaml.load_file(file.path(CONF_DIR, "seeds.yaml"))

THRESH  <- PARAMS$thresholds
PREP    <- PARAMS$preprocessing
FEATSEL <- PARAMS$feature_selection
MODEL   <- PARAMS$model
COMPUTE <- PARAMS$compute
COHORTS <- PARAMS$cohorts
TRAIN_COHORT      <- PARAMS$project$training_cohort
VALIDATION_COHORTS <- PARAMS$project$validation_cohorts

set.seed(SEED$global)

## ---- canonical directories (data exists once per file) --------------------
DIR_DATA      <- file.path(UC4_ROOT, "data")
DIR_RAW       <- file.path(DIR_DATA, "raw")
DIR_PROCESSED <- file.path(DIR_DATA, "processed")
DIR_DERIVED   <- file.path(DIR_DATA, "derived")
DIR_TABLES    <- file.path(DIR_DATA, "output")
DIR_METADATA  <- file.path(UC4_ROOT, "metadata")
DIR_FIGURES   <- file.path(UC4_ROOT, "figures")
DIR_MODEL     <- file.path(UC4_ROOT, "model")
DIR_LIB       <- file.path(UC4_ROOT, "lib")
DIR_SOURCE    <- file.path(UC4_ROOT, "figures", "source_data")

for (d in c(DIR_DATA, DIR_RAW, DIR_PROCESSED, DIR_DERIVED, DIR_TABLES,
            DIR_METADATA, DIR_FIGURES, DIR_MODEL, DIR_LIB))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- small utilities ------------------------------------------------------
msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
                        paste0(..., collapse = ""), "\n", sep = "")

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

cohort_dir <- function(cohort) ensure_dir(file.path(DIR_PROCESSED, cohort))

n_cores <- function() max(1L, min(COMPUTE$cores, parallel::detectCores() - 1L))

## ---- table I/O (all TSV, LF endings, row names preserved) -----------------
write_tsv <- function(x, path, row.names = FALSE, quote = FALSE, gz = FALSE) {
  ensure_dir(dirname(path))
  if (gz) {
    con <- gzfile(path, "w")
    on.exit(close(con))
    utils::write.table(x, con, sep = "\t", quote = quote,
                       row.names = row.names, col.names = TRUE)
  } else {
    utils::write.table(x, path, sep = "\t", quote = quote,
                       row.names = row.names, col.names = TRUE)
  }
  invisible(path)
}

read_expr <- function(path, gz = grepl("\\.gz$", path)) {
  dt <- data.table::fread(path, sep = "\t", header = TRUE,
                          data.table = FALSE, check.names = FALSE)
  rn <- dt[[1]]
  m <- as.matrix(dt[, -1, drop = FALSE])
  rownames(m) <- rn
  storage.mode(m) <- "double"
  m
}

read_group <- function(path) {
  g <- utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE,
                         colClasses = "character")
  if (ncol(g) >= 2) {
    colnames(g)[1:2] <- c("Sample", "Type")
    rownames(g) <- g$Sample
  }
  g
}

align_group <- function(expr, group) {
  common <- intersect(colnames(expr), group$Sample)
  expr <- expr[, common, drop = FALSE]
  group <- group[match(common, group$Sample), , drop = FALSE]
  list(expr = expr, group = group)
}

## ---- expression-matrix helpers --------------------------------------------
first_symbol <- function(x) {
  s <- strsplit(as.character(x), ";", fixed = TRUE)
  vapply(s, function(z) z[1], character(1))
}

read_annotation <- function(path, id_col = "GeneID", symbol_col = "Symbol") {
  ann <- data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE)
  ann <- ann[, c(id_col, symbol_col), drop = FALSE]
  colnames(ann) <- c("id", "symbol")
  ann$symbol <- first_symbol(ann$symbol)
  ann
}

map_probes_to_symbols <- function(expr, ann) {
  idx <- match(rownames(expr), ann$id)
  sym <- ann$symbol[idx]
  ok <- !is.na(sym) & nzchar(sym)
  expr <- expr[ok, , drop = FALSE]
  rownames(expr) <- toupper(sym[ok])
  expr
}

collapse_duplicate_symbols <- function(expr) {
  limma::avereps(expr, ID = rownames(expr))
}

normalise_expression <- function(expr, log2 = TRUE, impute = TRUE,
                                 pseudo = PREP$log2_pseudocount, k = PREP$impute_k) {
  if (log2) {
    expr[expr < 0] <- 0
    expr <- log2(expr + pseudo)
  }
  if (impute && anyNA(expr)) expr <- impute::impute.knn(expr)$data
  if (identical(PREP$normalization, "quantile"))
    expr <- limma::normalizeBetweenArrays(expr, method = "quantile")
  expr
}

qc_stats <- function(expr) {
  data.frame(sample = colnames(expr),
             min = apply(expr, 2, min), q25 = apply(expr, 2, quantile, 0.25),
             median = apply(expr, 2, median), q75 = apply(expr, 2, quantile, 0.75),
             max = apply(expr, 2, max),
             missing = colSums(is.na(expr)), row.names = NULL)
}

matrix_cols <- function(expr, genes) expr[intersect(genes, rownames(expr)), , drop = FALSE]

## ---- plotting / model helpers ---------------------------------------------
save_pdf <- function(plot, path, width = 8, height = 6) {
  ensure_dir(dirname(path))
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off())
  if (inherits(plot, "ggplot")) print(plot) else plot
  invisible(path)
}

predict_uc <- function(expr, fit = NULL) {
  genes <- MODEL$genes
  x <- matrix_cols(expr, genes)
  lp <- as.numeric(MODEL$coefficients["(Intercept)"]) +
        colSums(sweep(x, 1, vapply(genes, function(g) MODEL$coefficients[[g]], numeric(1)), "*"))
  1 / (1 + exp(-lp))
}

risk_score <- function(prob) ifelse(is.na(prob), NA_real_, as.numeric(prob >= MODEL$cutoff))

require_raw <- function(cohort, file) {
  p <- file.path(DIR_RAW, cohort, file)
  if (!file.exists(p)) stop("Missing raw file ", p,
                            ". Run R/00_download_data.R first (or set UC4_RAW).", call. = FALSE)
  p
}

## ---- source the remaining helper modules ----------------------------------
for (f in c("01_preprocess.R", "02_diff.R", "03_enrichment.R",
            "04_featuresel.R", "05_model.R", "06_immune.R")) {
  p <- file.path(UC4_ROOT, "R", f)
  if (file.exists(p)) source(p, local = FALSE)
}

msg("UC-4 root: ", UC4_ROOT)
msg("config: ", length(COHORTS), " cohorts, training = ", TRAIN_COHORT)
