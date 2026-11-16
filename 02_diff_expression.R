# =============================================================================
#  02_diff_expression.R
#  Main limma DEG analysis + two sensitivity analyses (duplicateCorrelation and
#  one-biopsy-per-patient) via the single run_deg_analysis() in R/02_diff.R,
#  followed by gene-set intersections (02.genesets.R) and GO/KEGG enrichment
#  (the clusterProfiler blocks of 03.venn.R / 04String.R).
#  Output: data/output/deg/, data/output/gene_sets/, data/output/enrichment/
# =============================================================================
source(file.path(Sys.getenv("UC4_HOME", getwd()), "utils.R"))
suppressPackageStartupMessages({ library(limma); library(data.table) })

expr <- read_expr(file.path(cohort_dir(TRAIN_COHORT), "expr_normalized.tsv"))
grp  <- read_group(file.path(cohort_dir(TRAIN_COHORT), "group.tsv"))
sp   <- read_group(file.path(cohort_dir(TRAIN_COHORT), "sample_patient_group.tsv"))
patients <- stats::setNames(sp$Patient, sp$Sample)

msg("running main + sensitivity DEG fits")
deg <- run_deg_analysis(expr, grp, patients = patients)

msg("gene-set intersections")
arg_deg <- gene_set_intersections()
enrichment(arg_deg, file.path(DIR_TABLES, "enrichment"), "ARG_DEG_intersection")
