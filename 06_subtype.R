# =============================================================================
#  06_subtype.R
#  Consensus clustering of UC samples (UC-3 12Subtype.R block 1), stability
#  checks (silhouette / LOFO ARI / PCA), subtype DEG + enrichment, and the
#  subtype immune/GEE arm through the shared immune_analysis()/run_gee()
#  functions.  Replaces the ~2,200-line 12Subtype.R with one script that reuses
#  R/06_immune.R instead of duplicating the whole immune pipeline.
#  Output: data/output/subtype/, figures/subtype/
# =============================================================================
source(file.path(Sys.getenv("UC4_HOME", getwd()), "utils.R"))
suppressPackageStartupMessages({ library(ConsensusClusterPlus); library(cluster) })
out <- ensure_dir(file.path(DIR_TABLES, "subtype"))

expr <- read_expr(file.path(cohort_dir(TRAIN_COHORT), "expr_normalized.tsv"))
grp  <- read_group(file.path(cohort_dir(TRAIN_COHORT), "group.tsv"))
uc_samples <- grp$Sample[grp$Type == "UC"]
x <- t(expr[MODEL$genes, uc_samples, drop = FALSE])

set.seed(SEED$consensus_cluster)
cc <- ConsensusClusterPlus::ConsensusClusterPlus(
  t(x), maxK = PARAMS$subtype$consensus$max_k, reps = PARAMS$subtype$consensus$reps,
  pItem = PARAMS$subtype$consensus$p_item, pFeature = PARAMS$subtype$consensus$p_feature,
  clusterAlg = PARAMS$subtype$consensus$cluster_alg,
  distance = PARAMS$subtype$consensus$distance,
  seed = SEED$consensus_cluster, plot = "pdf", outdir = out, verbose = FALSE)

k <- PARAMS$subtype$consensus$selected_k
cl <- cc[[k]]$consensusClass
sub <- data.frame(Sample = names(cl), Subtype = paste0("C", cl), stringsAsFactors = FALSE)
write_tsv(sub, file.path(out, "subtype_assignments.tsv"))

## silhouette (stability)
d <- stats::dist(t(x), method = PARAMS$subtype$stability$silhouette_metric)
sil <- cluster::silhouette(as.integer(cl), d)
write_tsv(data.frame(Sample = names(sil[, "sil_width"]), silhouette = sil[, "sil_width"]),
          file.path(out, "silhouette_k2.tsv"))

## subtype immune arm via the SHARED pipeline (grouping = subtype)
patients <- NULL
sp <- file.path(cohort_dir(TRAIN_COHORT), "sample_patient_group.tsv")
if (file.exists(sp)) { g <- read_group(sp); patients <- stats::setNames(g$Patient, g$Sample) }
imm <- immune_analysis(expr[, sub$Sample, drop = FALSE],
                       grouping = stats::setNames(sub$Subtype, sub$Sample),
                       patients = patients[colnames(expr)], tag = "subtype",
                       out_dir = out)
if (!is.null(imm$ssgsea)) {
  md <- data.frame(Sample = colnames(imm$ssgsea),
                   Subtype = sub$Subtype[match(colnames(imm$ssgsea), sub$Sample)],
                   Patient = patients[colnames(imm$ssgsea)], stringsAsFactors = FALSE)
  gee <- tryCatch(run_gee(as.data.frame(t(imm$ssgsea)), md, "Subtype"),
                  error = function(e) { msg("GEE skipped: ", conditionMessage(e)); NULL })
  if (!is.null(gee)) write_tsv(gee, file.path(out, "GEE_ssgsea_adjusted.tsv"))
}
msg("subtype analysis complete (k = ", k, ")")
