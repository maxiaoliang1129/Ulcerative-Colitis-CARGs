# Deduplication report

## Method

Every one of the 361 files in UC-3 was read and hashed with **MD5** of its full byte content. Files sharing an MD5 hash are byte-identical regardless of their name or folder. For each hash group one canonical copy was retained; the others are recorded as duplicates. No file was modified, and UC-3 was never changed.

## Summary

| Metric | Value |
|---|---:|
| Source files (UC-3) | 361 |
| Unique files (MD5) | 228 |
| Duplicate copies removed | 133 |
| Total size of all copies | 18.16 GB |
| Size retained (unique) | 4.59 GB |
| Space saved | 13.56 GB (74.7%) |
| Hash groups with duplicates | 55 |

## Largest duplicate groups

| Copies | Canonical file in UC-4 | Saved |
|---:|---|---:|
| 17 | data/processed/GSE193677/group_all_samples.txt | 0.4 MB |
| 16 | data/derived/model/final_model_genes.txt | 0.0 MB |
| 12 | data/processed/GSE193677/expr_normalized.txt | 10165.6 MB |
| 9 | data/derived/deg/diffGeneExp_p0.05.txt | 141.7 MB |
| 7 | data/derived/subtype/NMF_result__748b0451.csv | 0.1 MB |
| 5 | data/derived/ppi/mcode1.txt | 0.0 MB |
| 5 | data/derived/gene_sets/diff_intergene.txt | 0.0 MB |
| 5 | data/processed/GSE87466/expr_normalized.txt | 148.3 MB |
| 5 | data/processed/GSE128682/group.txt | 0.0 MB |
| 5 | data/processed/GSE128682/expr_normalized.txt | 42.7 MB |
| 5 | data/processed/GSE87466/group.txt | 0.0 MB |
| 4 | data/derived/ppi/Cluster1_node.csv | 0.0 MB |
| 3 | data/raw/GPL13158/GPL13158-5065.txt | 117.3 MB |
| 3 | data/processed/GSE193677/clinical_metadata.txt | 0.2 MB |
| 3 | data/derived/immune/uniq.symbol.txt | 1848.3 MB |
| 3 | data/derived/model/final_model_genes__66be95bf.txt | 0.0 MB |
| 3 | data/processed/GSE206285/group.txt | 0.0 MB |
| 3 | data/processed/GSE206285/expr_normalized.txt | 389.2 MB |
| 3 | data/processed/GSE193677/phenotype_selected.csv | 10.8 MB |
| 2 | data/raw/GSE128682/GSE128682_DESeq2_normalized_counts.tab | 4.1 MB |

## Notes

- The training `normalize.txt` (924 MB) appeared 12 times and is now one file at `data/processed/GSE193677/expr_normalized.txt`.
- Two files with the same name but different content (e.g. the 93-gene vs 34-gene `diff_intergene.txt`) were **not** merged; they were disambiguated by their canonical destinations and, where necessary, an MD5 suffix.
- Two zero-byte files (`04String/00string/Cluster1 (Score`, `11Immue/04cibersort/boxplot.pdf`) share the empty-file MD5 and were kept once.
- The full source-to-destination mapping is in `metadata/file_manifest.tsv`.
