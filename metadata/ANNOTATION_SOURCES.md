# Annotation sources — UC-4

| Resource | File | Origin | Used for |
|---|---|---|---|
| NCBI GRCh38.p13 gene annotation | `data/raw/GSE193677/Human.GRCh38.p13.annot.tsv` | GEO GSE193677 supplementary | GeneID → Symbol mapping for RNA-seq counts |
| GPL13158 platform annotation | `data/raw/GPL13158/GPL13158-5065.txt` | GEO GPL13158 | Probe → Gene symbol for GSE87466 / GSE206285 |
| GEO series matrices | `data/raw/<cohort>/*series_matrix*` | GEO (GSE193677/87466/206285) | Phenotype and group labels |
| DESeq2 normalised counts | `data/raw/GSE128682/GSE128682_DESeq2_normalized_counts.tab` | GEO GSE128682 supplementary | Expression for validation cohort |
| Autophagy gene set | `metadata/gene_sets/ARGs.txt` | original UC-3 (published list) | ARG ∩ DEG intersection |
| WGCNA hub genes | `metadata/gene_sets/wgcna_hubs.txt` | original UC-3 (fixed input) | candidate gene pool |
| Immune gene set | `metadata/gene_sets/immune.gmt` | original UC-3 | ssGSEA scoring |
| GO BP gene set | `metadata/gene_sets/c5.go.bp.v2024.1.Hs.symbols.gmt` | MSigDB v2024.1 | single-gene GSEA |
| CIBERSORT signature | `metadata/external/cibersort_signature.txt` | original UC-3 (CIBERSORT v1.03) | immune deconvolution |
| OrgDb | R package `org.Hs.eg.db` 3.19.1 | Bioconductor | Entrez / GO identifiers |

Gene-symbol mapping convention: the first symbol of a semicolon-separated field
is used (`R/00_config.R::first_symbol()`), duplicates are collapsed with
`limma::avereps()`. Symbols are upper-cased consistently across all cohorts.
