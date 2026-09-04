# SNE construction

This directory contains the analyses used to construct and inspect social niche
embeddings (SNEs) for the gut microbiome.

## 1. Scientific question

The analysis asks how co-occurrence choices affect the learned microbial representation and how the resulting SNE structure relates to phylogeny, taxonomy, traits, and prevalence. The code is organised into three parts:

1. construction of GloVe/SNE embeddings from microbial co-occurrence data;
2. comparison of co-occurrence metrics across study-level embeddings;
3. visualisation of SNE structure with phylogeny, taxonomy, traits, and prevalence.

The scripts use the `membed` package for co-occurrence calculation and GloVe
training.

## 2. Data

| Resource | Path | Used by |
|---|---|---|
| Full-cohort BIOM table | `../../data/gut_pretraining.biom` | `sne_training/` |
| Shared SNE embedding | `../../data/social_niche_embedding_100.txt` | `embedding_overview/` |
| Phylogenetic embedding | `../../data/phylo_embed_PCA_100.txt` | `SNE_t_SNE.ipynb` |
| Taxonomy and tree | `../../data/taxmap_slv_ssu_ref_nr_138.2.txt`, `../../data/SSURefNR99_1200_slv_138_2_subset.tre` | overview scripts |
| Study embeddings | `data/cooccurrence_metric_comparison/difference_study_training/embeding_list/` | metric comparison |
| Prevalence table | `data/embedding_overview/OTU_prevalence_abundance.csv` | overview scripts |

The checked-in full-cohort BIOM file is a Git LFS pointer. A real BIOM file is
needed to run the training scripts.

## 3. Scripts and notebooks

### SNE training

| File | Purpose |
|---|---|
| `sne_training/train_glove_from_source.sh` | Builds the feature dictionary, co-occurrence matrix, x-max file, and 100-dimensional GloVe embedding from a BIOM table. |
| `sne_training/train_glove_binary.sh` | Runs the corresponding binary-data training job. |
| `sne_training/train_glove_embedding_size.sh` | Repeats the training workflow with a selected embedding dimension. |

The wrappers write intermediates, embeddings, and logs to
`sne_training/results/`.

### Co-occurrence metric comparison

| File | Purpose |
|---|---|
| `cooccurrence_metric_comparison/train_metric_embeddings.sh` | Reads each BIOM path and metric from `param_matrix.txt` and trains one 100-dimensional embedding per study/metric pair. |
| `cooccurrence_metric_comparison/train_metric_embeddings_binary.sh` | Runs the same array using the binary-input parameter list in `param_matrix_binary.txt`. |
| `cooccurrence_metric_comparison/compute_similarity_matrices.py` | Converts each embedding file into an OTU-by-OTU cosine-similarity matrix. |
| `cooccurrence_metric_comparison/run_mantel_tests.py` | Groups matrices by metric, compares every study matrix with the metric baseline, and writes Mantel R and permutation P-values. |
| `cooccurrence_metric_comparison/plot_mantel_heatmap.ipynb` | Converts the Mantel result table into the metric-by-study heatmap. |
| `cooccurrence_metric_comparison/plot_metric_networks.R` | Samples 20 OTUs from one AGP sample, calculates eight metric-specific edge weights, and displays the network comparison. |

### Embedding overview

| File | Purpose |
|---|---|
| `embedding_overview/SNE_t_SNE.ipynb` | Computes cosine-distance t-SNE coordinates for the SNE and phylogenetic embeddings, including the BugBase subset, and writes the coordinate tables used by the plots. |
| `embedding_overview/SNE_tree_overview.R` | Orders SNE dimensions, displays them alongside the phylogenetic tree, and saves the combined figure. |
| `embedding_overview/plot_tax_bugbase.R` | Colours SNE and phylogenetic t-SNE coordinates by taxonomy and selected BugBase traits. |
| `embedding_overview/plot_tsne_prevalence.py` | Colours a two-dimensional SNE projection by OTU prevalence. |

## 4. Main outputs

- `sne_training/results/` — training intermediates, embeddings, and logs;
- `cooccurrence_metric_comparison/results/similarity_matrix/` — cosine-similarity matrices;
- `cooccurrence_metric_comparison/results/mantel_test_results.csv` — Mantel-test results;
- `embedding_overview/results/tables/` — t-SNE coordinate tables;
- `embedding_overview/results/figures/` — overview figures.
