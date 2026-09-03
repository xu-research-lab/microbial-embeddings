# SNE construction

This analysis has three parts:

1. compare microbial co-occurrence metrics across cohorts;
2. train the final Social Niche Embedding (SNE);
3. compare SNE with phylogeny, prevalence, taxonomy, and microbial traits.

Run commands from the repository root. Python commands use the
`microbiome_deep` environment.

## Data

Shared inputs are stored in the repository-level `data/` directory:

| File | Use |
| --- | --- |
| `data/gut_pretraining.biom` | SNE pre-training table (14,093 features, 202,558 samples) |
| `data/table_gut_all.biom` | Complete gut table before disease-sample removal |
| `data/metadata_gut_all.tsv` | Sample metadata for the complete gut table |
| `data/social_niche_embedding_100.txt` | Final 100-dimensional SNE |
| `data/Embedding_list/*.txt` | Full-cohort embeddings for the eight metrics |
| `data/phylo_embed_PCA_100.txt` | 100-dimensional phylogenetic embedding |
| `data/taxmap_slv_ssu_ref_nr_138.2.txt` | SILVA taxonomy |
| `data/SSURefNR99_1200_slv_138_2_subset.tre` | Phylogenetic tree |

Analysis-specific cohort tables and embedding runs remain under
`analysis/sne_construction/data/`. New computed outputs go under the
corresponding `results/` directory.

## Train an embedding

The shared Slurm entry point runs feature counting, co-occurrence calculation,
`x_max` calculation, and GloVe training:

```bash
sbatch analysis/sne_construction/sne_training/run_glove.sh
```

With no arguments it uses the adopted settings documented by the project:
`gut_pretraining.biom`, `abundance_percentile`, `x_max=80`, 100 dimensions,
learning rate 0.05, and 100 iterations. Results are written to
`analysis/sne_construction/sne_training/results/abundance_percentile_100/`.
The script never replaces `data/social_niche_embedding_100.txt`.

Alternative runs must state their changed parameter and output directory:

```bash
sbatch analysis/sne_construction/sne_training/run_glove.sh \
  data/gut_pretraining.biom jaccard 50 \
  analysis/sne_construction/sne_training/results/jaccard_50 28
```

Arguments are `BIOM METRIC EMBEDDING_SIZE OUTPUT_DIR CPUS COOCCURRENCE_DIR
XMAX_PERCENTILE`. The last two are optional; provide a directory containing
`table.co` and `feature-dict.txt` or `feature-dict.csv` to reuse an existing
co-occurrence calculation.

## Compare cohort embeddings

The retained study tables and embedding families are in
`analysis/sne_construction/data/cooccurrence_metric_comparison/`.
Use `run_glove.sh` for any missing cohort/metric run, then compute similarity
matrices and Mantel statistics:

```bash
conda activate microbiome_deep
python analysis/sne_construction/cooccurrence_metric_comparison/compute_similarity.py
python analysis/sne_construction/cooccurrence_metric_comparison/run_mantel_test.py
```

Default outputs:

- `cooccurrence_metric_comparison/results/similarity_matrices/`
- `cooccurrence_metric_comparison/results/mantel_test_results.csv`
- `cooccurrence_metric_comparison/results/figures/`

`compute_similarity.py` creates dense matrices and defaults to one worker to
avoid multiplying its memory demand. The non-rarefied cohort tables can be
regenerated from `data/table_gut_all.biom`. The retained rarefied cohort tables
remain under `data/cooccurrence_metric_comparison/difference_study_training/data_rarefaction/`; their source rarefied compendium is not in the repository.

## Embedding overview

`embedding_overview/SNE_t_SNE.ipynb` compares SNE with the phylogenetic
embedding and writes four retained coordinate tables to
`embedding_overview/results/tables/` for taxonomy and trait plots.
`plot_tsne_prevalence.py` writes
`embedding_overview/results/figures/tsne_prevalence.png`.
The R scripts write their figures to the same results tree.
