# Synthetic validation

This module generates consumer-resource-model (CRM) communities and evaluates
how well synthetic microbial embeddings recover the simulated niche structure.

Start in `analysis/synthetic_validation`. The visualization script is run
from `embedding_validation` because it reads and writes that directory's
`results/plot_ratios/` tree.

## Data

- `data/E_global.csv`: species-by-resource interaction matrix.
- `data/metrics_embeddings/`: five embedding replicates for each co-occurrence
  metric. The `rep1`-`rep5` suffixes are experimental replicates.
- `data/datasize_embeddings/`: five 100-dimensional embedding replicates for
  each training-set size.
- `data/datasize_bioms/datasize_bioms/`: five tracked BIOM subsets for each
  training-set size. Moving the local HDF5 copies is deferred and is not part
  of this pull request.
- `embedding_validation/results/plot_ratios/plot_csvs/`: retained ratio tables.
  Filename suffixes encode the high and low similarity thresholds, not versions.

## Analyses

```bash
cd analysis/synthetic_validation

Rscript crm_generation/run_CRM.r

mkdir -p results
conda run -n microbiome_deep python data_size_ablation/datasize.py

cd embedding_dimension_ablation
conda run -n microbiome_deep python embedding_dim_iter.py
cd ..

jupyter lab embedding_validation/

cd embedding_validation
Rscript ../visualization/box_plot_log.R
```

`run_CRM.r` writes `data/E_global.csv` and
`data/final_abundance_table.csv`. Data-size summaries and figures are written to
`results/`; embedding-dimension outputs are written to
`embedding_dimension_ablation/results/`; ratio plots are written to
`embedding_validation/results/plot_ratios/`.

The retained repository does not include the embedding-dimension training logs
or the raw `table.co` and `feature-dict.csv` files required by
`analysis_embedding.ipynb`. The notebook expects the latter two files under
`data/cooccurrence/`.
