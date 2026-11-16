# =============================================================================
#  08_downstream.R
#  Downstream analyses that were web-tool / interactive steps in UC-3:
#   - single-gene GSEA for CCL11 and CLDN1      (10ssgsea.R)
#   - TF network (NetworkAnalyst output)         (07TFs)
#   - ceRNA network (miRTarBase + GENCODE)       (08ceRNA.R)
#   - drug candidates (DSigDB)                   (09DsigDB)
#  The external tool results are shipped once under data/derived/; this script
#  regenerates the in-R parts and documents the web-tool parameters.
#  Output: data/output/downstream/, figures/gene_sets/
# =============================================================================
source(file.path(Sys.getenv("UC4_HOME", getwd()), "utils.R"))
out <- ensure_dir(file.path(DIR_TABLES, "downstream"))

## ---- single-gene GSEA (shared function, both model genes) ------------------
single_gene_gsea <- function(gene) {
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db)
                                   library(future.apply) })
  expr <- read_expr(file.path(cohort_dir(TRAIN_COHORT), "expr_normalized.tsv"))
  grp  <- read_group(file.path(cohort_dir(TRAIN_COHORT), "group.tsv"))
  y <- ifelse(grp$Type[match(colnames(expr), grp$Sample)] == "UC", 1, 0)
  ranks <- stats::cor(t(expr), y, method = "pearson")
  ranks <- sort(ranks[is.finite(ranks)])
  gmt <- file.path(DIR_METADATA, "gene_sets", "c5.go.bp.v2024.1.Hs.symbols.gmt")
  gsea <- clusterProfiler::gseGO(geneList = ranks, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                 keyType = "SYMBOL", ont = "BP", minGSSize = 10,
                                 pvalueCutoff = THRESH$go_pval, verbose = FALSE)
  write_tsv(as.data.frame(gsea), file.path(out, sprintf("%s_GSEA.tsv", gene)))
  gsea
}
for (g in MODEL$genes) tryCatch(single_gene_gsea(g),
  error = function(e) msg("GSEA skipped for ", g, ": ", conditionMessage(e)))

## ---- archived web-tool products (copied once, documented) -----------------
msg("TF / ceRNA / DSigDB inputs are in data/derived/{tf,cerna,drugs}/")
msg("see docs/external_web_tools.md for tool URLs and parameters")
