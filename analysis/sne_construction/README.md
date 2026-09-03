# SNE construction

This module groups co-occurrence metric comparison, SNE training, and embedding
overview analyses. Reusable implementation is provided by the repository's
`membed/` package; the files here are retained analysis wrappers and notebooks.

Python, R, and notebook paths are relative to `analysis/sne_construction/`; run
those analyses from that directory. The legacy Slurm wrappers retain their original
working-directory assumptions because their parameter files are unavailable.

| Subanalysis | Core files | Inputs | Outputs | Status |
| --- | --- | --- | --- | --- |
| Co-occurrence metric comparison | `cooccurrence_metric_comparison/glove_train_source.sh`, `glove_train_source_binary.sh`, `compute_similarity.py`, `run_metal_test.py`, `analysis_metric.ipynb`, `Comparative_analysis.R` | BIOM subsets and metric embedding families under `data/cooccurrence_metric_comparison/` | retained heatmap under `cooccurrence_metric_comparison/results/figures/`; Mantel and similarity families require path/provenance reconciliation | Legacy source; adopted metric/run unknown |
| SNE training | `sne_training/glove_trainning_source.sh`, `glove_trainning_with_embedding_size.sh`, `glove_trianning_binary.sh` | `../../data/gut_pretraining.biom`, parameter matrices, and optional co-occurrence intermediates | logs and run-specific summaries belong under `sne_training/results/`; the retained shared embedding is `../../data/social_niche_embedding_100.txt` | Legacy source; adopted invocation unknown |
| Embedding overview | `embedding_overview/SNE_t_SNE.ipynb`, `SNE_tree_overview.R`, `plot_tax_bugbase.R`, `plot_tsne_prevalence.py` | shared SNE, tree, and taxonomy under `../../data/`; prevalence under `data/embedding_overview/`; traits under `../traits/data/`; t-SNE tables under `embedding_overview/results/tables/` | figures belong under `embedding_overview/results/figures/` | Legacy source; adopted figures unknown |

`../../data/social_niche_embedding_100.txt` is a repository-level shared artifact
used by this analysis and downstream analyses. The natural dependency is compendium data -> metric/training artifacts -> shared
SNE -> overview and downstream analyses.

The retained wrappers do not identify the adopted parameter matrix, metric,
`xmax`, seed, repetition, or run manifest. Several source paths remain
external or absent; reorganization alone does not establish a reproducible
training run.
