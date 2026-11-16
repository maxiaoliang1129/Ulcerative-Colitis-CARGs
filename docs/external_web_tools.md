# External web tools and archived analyses

The following UC-3 steps were performed with external web tools or
interactively and are shipped as their result files (one canonical copy each).
Their in-R preparation/regeneration is in `08_downstream.R`.

| Tool / analysis | UC-3 folder | Result file in UC-4 | Notes |
|---|---|---|---|
| STRING protein interactions | `04String/00string/` | `data/derived/ppi/string_interactions_short.tsv` | Homo sapiens; confidence score from `config/params.yaml` (`mcode_string_score: 0.4`) |
| Cytoscape MCODE | `04String/00string/` | `data/derived/ppi/mcode1.txt`, `Cluster*_node.csv` | MCODE version/parameters were not recorded in UC-3 — to be stated in Methods |
| NetworkAnalyst TF network | `07TFs/` | `data/derived/tf/networkanalyst_0.json`, `NetworkAnalyst.json.{sif,cyjs,svg}` | TF–gene interactions for the two model genes |
| miRTarBase (hsa_MTI) + GENCODE v46 | `08ceRNA/` | `metadata/external/hsa_MTI.xlsx`, `data/derived/cerna/ceRNA.csv` | lncRNA–miRNA–mRNA network intersected with CCL11/CLDN1 |
| DSigDB (Enrichr) | `09DsigDB/` | `data/derived/drugs/DSigDB_table_select.xlsx` | candidate compounds for the model genes |
| CIBERSORT (v1.03, vendored) | `11Immue/04cibersort/` | `lib/CIBERSORT.R`, `metadata/external/cibersort_signature.txt` | 100 permutations, SVR-based deconvolution |
| WGCNA hub genes | `02genesets/ARGS/` | `metadata/gene_sets/wgcna_hubs.txt` | fixed input; no WGCNA script existed in UC-3 |

All network graph files (`.sif`, `.cyjs`, `.cys`, `.svg`) are single copies under
`data/derived/ppi/`, `data/derived/tf/` or `data/derived/cerna/`.
