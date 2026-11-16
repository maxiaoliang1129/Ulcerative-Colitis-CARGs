# =============================================================================
#  00_download_data.R
#  Fetch the four public GEO series + the shared GPL13158 platform table into
#  data/raw/.  Files already present are skipped.  Run once before 01.
#  Replaces the download portions of the four data*.R scripts in UC-3.
# =============================================================================
source(file.path(Sys.getenv("UC4_HOME", getwd()), "utils.R"))
suppressPackageStartupMessages(library(GEOquery))

geo_base <- "https://ftp.ncbi.nlm.nih.gov/geo/series"
files <- c(
  "GSE193677" = c("GSE193677_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                  "GSE193677_family.soft.gz", "GSE193677_series_matrix.txt.gz",
                  "Human.GRCh38.p13.annot.tsv.gz"),
  "GSE128682" = c("GSE128682_DESeq2_normalized_counts.tab.gz"),
  "GSE87466"  = c("GSE87466_series_matrix.txt.gz"),
  "GSE206285" = c("GSE206285_series_matrix.txt.gz")
)

for (cohort in names(files)) {
  outdir <- ensure_dir(file.path(DIR_RAW, cohort))
  series <- substr(cohort, 1, 6)             # GSE1936 -> GSE1936nnn grouping
  base <- file.path(geo_base, paste0(substr(cohort, 1, 6), "nnn"), cohort, "suppl")
  for (f in files[[cohort]]) {
    dest <- file.path(outdir, f)
    if (file.exists(dest)) { msg("skip (exists): ", f); next }
    msg("downloading ", cohort, "/", f)
    tryCatch(utils::download.file(file.path(base, f), dest, mode = "wb", quiet = TRUE),
             error = function(e) msg("  FAILED: ", conditionMessage(e)))
    if (grepl("\\.gz$", f) && file.exists(dest) && !file.exists(sub("\\.gz$", "", dest)))
      R.utils::gunzip(dest, remove = FALSE)
  }
}
# shared platform annotation, stored once
gpldir <- ensure_dir(file.path(DIR_RAW, "GPL13158"))
if (!file.exists(file.path(gpldir, "GPL13158-5065.txt")))
  msg("place GPL13158-5065.txt in data/raw/GPL13158/ (from GEO GPL13158)")
msg("done. Set UC4_RAW to point at an existing raw-data directory if downloads fail.")
