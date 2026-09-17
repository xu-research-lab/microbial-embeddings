# The `data/` directory

Every input file the analyses read lives here. Nothing in this folder is
generated on the fly: each file was either downloaded from a public database or
produced by the pipelines in [`analysis/resources/`](../analysis/resources/),
which take weeks of compute.

**First check:** these files are stored with [Git LFS](https://git-lfs.com). If
`gut_pretraining.biom` is only ~134 bytes, LFS did not run — go back to step 0
of the [main README](../README.md).

## Which file do I need?

| I want to | Use |
|---|---|
| use the SNE vectors from the paper | `social_niche_embedding_100.txt` |
| train my own SNEs from scratch | `gut_pretraining.biom` |
| run the disease classifier | `social_niche_embedding_removing_disease_samples_100.txt` + `metadata_disease_classification.tsv` |
| compare SNEs against other embeddings | `phylo_embed_PCA_100.txt`, `dnabert2_16s_embedding_reduced_100.txt` |
| look up what an OTU is | `taxmap_slv_ssu_ref_nr_138.2.txt` |

## The two formats you will meet

### 1. Embedding files (`*.txt`)

Plain text, one taxon per row, values separated by a single space, **no header**:

```
AAAA02020714.1.1202 -0.309345 -0.002363 -0.273315 ...   (100 numbers)
```

The first column is the OTU ID — a SILVA accession such as
`AAAA02020714.1.1202`. The remaining columns are its vector. The BIOM tables use
the same IDs as feature IDs, so tables and embeddings join directly.

Files trained with GloVe carry one extra row, `<unk>`. It is the fallback vector
for a taxon the model never saw; drop it if you only want real taxa.

```python
import pandas as pd
emb = pd.read_csv("social_niche_embedding_100.txt", sep=" ", header=None, index_col=0)
emb = emb.drop("<unk>", errors="ignore")   # 14,093 taxa x 100 dimensions
```

### 2. Abundance tables (`*.biom`)

[BIOM](https://biom-format.org) HDF5 files: rows are taxa, columns are samples,
values are read counts.

```python
import biom
table = biom.load_table("gut_pretraining.biom")
print(table.shape)          # (taxa, samples)
df = table.to_dataframe(dense=True)
```

## File list

### Abundance tables

| File | Size | What it is |
|---|---|---|
| `table_gut_all.biom` | 173M | The full compendium: **14,093 taxa x 210,090 samples**, merged from 643 public 16S studies and clustered against SILVA 138.2 NR99 at 97% identity. Everything else starts from this file. |
| `gut_pretraining.biom` | 164M | `table_gut_all.biom` with every sample of the disease benchmark removed: **14,093 taxa x 202,558 samples**. Train SNEs on this one, so the embedding never sees the samples it is later tested on. Built by `analysis/resources/16S_database/get_pretraining_datasets.ipynb`. |

### Sample metadata

| File | Size | What it is |
|---|---|---|
| `metadata_gut_all.tsv` | 17M | One row per sample of `table_gut_all.biom` (210,090 rows). Columns: `sample_id`, `project`, `instrument`, `geo_loc_name`, `region`, `seq_region` (the 16S variable region), `sample_type`. |
| `metadata_gut_all.txt` | 17M | Byte-identical copy of the file above, kept for older scripts. |
| `metadata_disease_classification.tsv` | 1.6M | The 10,276 samples of the disease benchmark, from 54 case-control studies covering 13 diseases (IBD, PD, OB, T2DM, IBS, CRC, SZ, MS, BD, ASD, GD, AS, CAD). Main columns: `sample` (the sample ID), `study`, `group` (the 0/1 label the classifier trains on), `diagnosis`, `disease_name_ab`, plus host and sequencing details. The ID column is named `sample`, not `sample_id`, so pass `--sample-id-col sample`. |

### Embeddings from the paper

| File | Size | What it is |
|---|---|---|
| `social_niche_embedding_100.txt` | 13M | **The main result.** 100-dimensional SNEs for all 14,093 taxa, trained on the full compendium with the `abundance_percentile` metric. Use it for any descriptive analysis. |
| `social_niche_embedding_removing_disease_samples_100.txt` | 13M | The same, but trained on `gut_pretraining.biom`. **Use this one for disease classification** — it is what keeps the benchmark clean. |
| `social_niche_embedding_removing_disease_samples_100_shuffled.txt` | 13M | The file above with each of the 100 dimensions permuted independently across taxa (`analysis/Disease_classification_loo/SNEs_shuffle.py`, seed 5), so every dimension keeps its values but the taxon-level structure is destroyed. A negative control: a model given these should lose its advantage. |

### Baseline embeddings

Alternative ways to represent the same taxa, used to show what the SNEs add.

| File | Size | What it is |
|---|---|---|
| `phylo_embed_PCA_100.txt` | 28M | Phylogeny instead of ecology: pairwise distances on the SILVA tree, reduced to 100 PCA dimensions. |
| `dnabert2_16s_embedding.txt` | 137M | Raw DNABERT-2 output, 768 dimensions per 16S sequence, computed from the sequence alone. |
| `dnabert2_16s_embedding_reduced_100.txt` | 27M | The file above reduced to 100 PCA dimensions, so it can be dropped in wherever the SNEs are used. |
| `other_embedding.ipynb` | 24K | The notebook that produced the three files above. |

### `Embedding_list/` — one SNE per co-occurrence metric

Eight 100-dimensional embeddings trained on the same data, each with a different
definition of "these two microbes co-occur". They back the metric comparison in
the paper.

| File | Metric |
|---|---|
| `abundance_percentile_100.txt` | abundance-weighted, percentile-ranked — **the metric chosen for the paper** (identical to `social_niche_embedding_100.txt`) |
| `abundance_totalsum_100.txt` | abundance-weighted, total-sum normalised |
| `braycurtis_percentile_100.txt`, `braycurtis_totalsum_100.txt` | Bray-Curtis similarity |
| `russell_rao_weight_100.txt` | rank-weighted co-occurrence |
| `russell_rao_100.txt`, `jaccard_100.txt`, `faith_100.txt` | presence/absence only |

The three presence/absence files cover 12,606 taxa instead of 14,093; the
remaining taxa get no vector under those metrics.

### SILVA 138.2 reference files

| File | Size | What it is |
|---|---|---|
| `taxmap_slv_ssu_ref_nr_138.2.txt` | 67M | Tab-separated SILVA taxonomy, 510,495 rows. Columns: `primaryAccession`, `start`, `stop`, `path` (the full `Bacteria;Bacillota;...` lineage), `organism_name`, `taxid`. To name an OTU, match the accession part of its ID — for `AB002518.1.1416`, look up `AB002518`. |
| `SSURefNR99_1200_slv_138_2_subset.tre` | 516K | Newick tree of the 14,093 OTUs, cut out of the SILVA reference tree. Source of every phylogenetic distance in the paper. |

### `projects/` — the raw material

1.5 GB, 643 subfolders, one per study, named by its archive accession (e.g.
`PRJEB13679/`). Each holds a single `DADA2.biom`: that study's ASV table as it
came out of DADA2, before SILVA clustering. Feature IDs are the ASV sequences
themselves and sample IDs are run accessions.

You need these only to rebuild the compendium from scratch. For everything else,
use `table_gut_all.biom`.

`projects/DADA2_table_summary.txt` lists all 643 studies with their ASV count,
sample count, and table density — a quick way to see what went in.

---

How each file was produced, and how to regenerate it, is documented in
[`analysis/resources/README.md`](../analysis/resources/README.md).
