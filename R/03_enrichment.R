# =============================================================================
#  R/03_enrichment.R  --  gene-set intersections & enrichment (shared)
#
#  One enrichment() function replaces the four copy-pasted GO/KEGG blocks in
#  UC-3/03.venn.R and UC-3/04String.R (and the 12Subtype GO/KEGG block).
#  KEGG is wrapped in tryCatch so an offline run still produces GO results.
# =============================================================================

#' GO / KEGG over-representation for one gene list.
#' @param genes   vector of gene symbols
#' @param out_dir where the four result tables are written
#' @param prefix  tag for the output file names
enrichment <- function(genes, out_dir, prefix) {
  ensure_dir(out_dir)
  genes <- unique(toupper(genes))
  if (!length(genes)) { msg("enrichment: empty gene list for ", prefix); return(invisible(NULL)) }
  suppressPackageStartupMessages({
    library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)
  })
  entrez <- as.character(org.Hs.eg.db::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                             keys = genes, column = "ENTREZID",
                                             keytype = "SYMBOL", multiVals = "first"))
  entrez <- entrez[!is.na(entrez)]

  go <- clusterProfiler::enrichGO(gene = entrez, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                  keyType = "ENTREZID", ont = "ALL",
                                  pvalueCutoff = THRESH$go_pval, qvalueCutoff = THRESH$go_qval,
                                  readable = TRUE)
  write_tsv(as.data.frame(go), file.path(out_dir, sprintf("%s_GO.tsv", prefix)))
  p <- clusterProfiler::dotplot(go, showCategory = 15) + ggplot2::theme_bw()
  save_pdf(p, file.path(out_dir, sprintf("%s_GO_dotplot.pdf", prefix)))

  kegg <- tryCatch({
      k <- clusterProfiler::enrichKEGG(gene = entrez, organism = "hsa",
                                       pvalueCutoff = THRESH$kegg_pval)
      write_tsv(as.data.frame(k), file.path(out_dir, sprintf("%s_KEGG.tsv", prefix)))
      save_pdf(clusterProfiler::dotplot(k, showCategory = 15) + ggplot2::theme_bw(),
               file.path(out_dir, sprintf("%s_KEGG_dotplot.pdf", prefix)))
      k
    }, error = function(e) {
      msg("KEGG not run for ", prefix, ": ", conditionMessage(e))
      writeLines(conditionMessage(e), file.path(out_dir, sprintf("%s_KEGG_NOT_RUN.txt", prefix)))
      NULL
    })
  list(go = go, kegg = kegg)
}

## ---- intersections between the published gene lists -----------------------
gene_set_intersections <- function(out_dir = file.path(DIR_TABLES, "gene_sets")) {
  ensure_dir(out_dir)
  read_list <- function(f) trimws(readLines(f, warn = FALSE))
  args  <- read_list(file.path(DIR_METADATA, "gene_sets", "ARGs.txt"))
  degs  <- read_list(file.path(DIR_METADATA, "gene_sets", "DEGs.txt"))
  wgcna <- read_list(file.path(DIR_METADATA, "gene_sets", "wgcna_hubs.txt"))
  arg_deg <- intersect(args, degs)
  writeLines(arg_deg, file.path(out_dir, "ARG_DEG_intersection.txt"))
  writeLines(intersect(arg_deg, wgcna), file.path(out_dir, "ARG_DEG_WGCNA_intersection.txt"))
  arg_deg
}

plot_intersections <- function(lists, path) {
  p <- UpSetR::upset(UpSetR::fromList(lists), order.by = "freq",Sets = c("ARGs","DEGs","WGCNA"))
  grDevices::pdf(path, width = 8, height = 6); on.exit(grDevices::dev.off()); print(p)
}
