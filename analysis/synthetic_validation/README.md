# Synthetic validation

This module contains CRM synthetic-community generation and validation analyses
associated with Fig. 1B-C, Extended Data Fig. 3A-B, and Supplementary Fig.
S3G-H. File-level figure adoption is not established.

| Subanalysis | Core files | Inputs | Outputs | Status |
| --- | --- | --- | --- | --- |
| CRM generation | `crm_generation/run_CRM.r` | embedded CRM parameters and package environment | shared `data/E_global_v46.csv`; expected final abundance table is absent | Candidate generator; adopted run unknown |
| Embedding validation | `embedding_validation/analysis_embedding.ipynb`, `analysis_distance_metrics.ipynb` | CRM data and metric embedding family under `data/metrics_embeddings/` | retained ratio CSVs under `embedding_validation/results/plot_ratios/plot_csvs/` | Exploratory notebooks; input/run provenance incomplete |
| Data-size ablation | `data_size_ablation/datasize.py` | retained BIOM and embedding families under `data/datasize_bioms/` and `data/datasize_embeddings/` | summaries and figures belong under `data_size_ablation/results/` | External input paths remain in code |
| Embedding-dimension ablation | `embedding_dimension_ablation/embedding_dim_iter.py` | retained logs under `embedding_dimension_ablation/results/embedding_dim_log/` | `plot_dimension.csv` and figures belong under the same result directory | External path remains in code |
| Visualization | `visualization/box_plot_log.R` | ratio CSVs produced by embedding validation | visualization outputs belong under `visualization/results/plot_ratios/` | Consumer path still requires reconciliation |

`data/E_global_v46.csv` is a shared CRM product consumed by multiple synthetic
subanalyses, so it remains at module level. Task-specific results remain with
the producing subanalysis. In particular, ratio CSVs belong to embedding
validation; the visualization task consumes them rather than copying them.

The final CRM abundance table, synthetic SNE training workflow, source tables
for the reported Mantel statistics, and adopted-run manifests are still absent.
