# Trait Analysis

This directory contains all phenotype-related analyses of the
[SNEs](https://github.com/xu-research-lab/microbial-embeddings) project. We ask
whether microbial **Social Niche Embeddings (SNEs)** — 100-dimensional vectors
pretrained from co-occurrence patterns in >210,000 human gut microbiomes —
encode the *phenotypic traits* of microbes (oxygen preference, Gram status,
cell shape, sporulation, motility, and metabolic capabilities), and whether
this trait structure is anything more than phylogenetic relatedness.

The central document of this analysis is
**[`traits_results.ipynb`](traits_results.ipynb)** (an R notebook). It reads
the results produced by the scripts below and assembles **Figure 2** and
**Extended Data Fig. 5** of the study (written to `results/`). This README
walks through the notebook panel by panel and explains what each script does,
which files it reads, and what it writes.

---

## 1. Scientific questions

1. **Analogy (Fig. 2A)** — Does the embedding satisfy the vector relation
   "facultatively anaerobic − aerobic ≈ anaerobic" (the same kind of arithmetic
   as "king − man + woman ≈ queen")?
2. **Trait structure beyond phylogeny (Fig. 2B, Ext. Fig. 5)** — Do taxa with
   the same trait cluster in SNE space? PLS-DA separates trait classes, and
   significance is judged against a *phylogenetic* null model (999 Brownian-motion
   threshold simulations on the SILVA tree), so a significant trait must be
   better separated than a random trait with the same degree of phylogenetic
   clustering.
3. **Generalization (Fig. 2C)** — Can a random forest read a trait off the
   embedding when trained on one database and tested on another, and when the
   test phylum was held out of training?
4. **Scaling (Fig. 2D)** — Does trait information grow with the size of the
   pretraining corpus?
5. **Embedding comparison (Fig. 2C)** — How do SNEs compare with DNABERT2 and
   phylogenetic-PCA embeddings for all of the above?

## 2. The trait data sources

Phenotypes come from **three independent databases**, so results do not depend
on a single annotation pipeline:

| Source | Kind of annotation | Raw output | Curated table used downstream |
|---|---|---|---|
| BugBase | Precalculated traits mapped from 16S sequences | `data/traits_bugbase.csv` |
| Traitar | Traits predicted from genomes | `data/trait_predcit.csv` |
| BacDive | Curated database records | — | `data/bacDive.csv` |

Cellular traits (Gram status, oxygen preference, cell shape, spore formation,
motility) are available in all three sources; metabolic traits (sugar
utilization: "Growth: Sugar" in Traitar, "assimilation" / "builds acid from"
in BacDive) come from Traitar and BacDive.

## 3. Shared input data

| Resource | Path | Description |
|---|---|---|
| SNE embedding | `../../data/social_niche_embedding_100.txt` | 100-d SNE vectors, one row per taxon |
| DNABERT2 embedding | `../../data/dnabert2_16s_embedding_reduced_100.txt` | DNABERT2 16S embedding, reduced to 100 d |
| Phylogenetic embedding | `../../data/phylo_embed_PCA_100.txt` | Phylogenetic embedding, PCA-reduced to 100 d |
| Subsample SNEs | `data/embedding/subset_table_{1w..16w}_100_{1..5}.txt` | SNEs retrained on 10k–160k-sample pretraining subsets, 5 replicates each |
| Taxonomy | `../../data/taxmap_slv_ssu_ref_nr_138.2.txt` | SILVA 138.2 taxonomy (used to map taxa to phyla) |
| Phylogenetic tree | `../../data/SSURefNR99_1200_slv_138_2_subset.tre` | SILVA tree subset (phylogenetic null model, phylolm) |

## 4. Scripts

### 4.1 Trait annotation (`traits_annotation/`)

| Script | Purpose | Reads | Writes |
|---|---|---|---|
| [`run_vsearch.sh`](traits_annotation/run_vsearch.sh) | Global-aligns the 16S query sequences against the BugBase reference database (Greengenes 99% OTU representatives) with `vsearch --usearch_global --id 0.99`. | `../data/feces_seq_16S_silva.fasta` (queries), `../data/99_otus.fasta` (database) | `../data/vsearch.out` (blast6out format) |
| [`Bugbase_predict.py`](traits_annotation/Bugbase_predict.py) | Transfers BugBase phenotype predictions onto the query sequences: keeps only vsearch hits that are **exact full-length matches** (query start = 1, 100% coverage, 100% identity) and relabels the BugBase trait table with the query IDs. Unmatched queries are dropped. | `../data/vsearch.out`, `Bugbase_database.txt` (BugBase trait table indexed by reference sequence ID), query FASTA for sequence lengths | `traits_predict_Bugbase.txt` (copy: `../data/traits_precalculated.txt`) |
| [`traitar_predict.sh`](traits_annotation/traitar_predict.sh) | SLURM job running `traitar phenotype` to predict phenotypes from genomes (Pfam 33.1 database, `from_nucleotides` mode, 28 CPUs). | `samples.txt` (genome list), Pfam database, genome files | `traits_predict_Traitar.csv` |


### 4.2 Trait–embedding analyses

| Script | Purpose | Reads | Writes |
|---|---|---|---|
| [`run_plsda.R`](run_plsda.R) | **PLS-DA + phylogenetic permutation test.** For every trait of every source, fits a PLS-DA (`ropls`, 2 predictive components, 5-fold CV) on the SNE embedding and tests the observed Q2/R2Y against 999 null traits simulated by a **Brownian-motion threshold model** on the SILVA tree (null traits have the same class counts and the same phylogenetic clustering as the observed trait). Stores the observed scores and statistics per trait. | `../../data/social_niche_embedding_100.txt`, `../../data/SSURefNR99_1200_slv_138_2_subset.tre`, `data/traits_bugbase.csv`, `data/trait_predcit.csv`, `data/bacDive.csv`, `data/agg_bac.csv` | `data/{bugbase,traitar,bacdive}/phyloperm_<trait>.rds` (scores, Q2, R2Y, pQ2/pR2Y, z-scores, class counts), `data/logs/error_*.log`, optionally `data/fig2b_phylo_perm.csv` |
| [`traits_predict.ipynb`](traits_predict.ipynb) | **Cross-database, cross-phylum trait prediction.** Trains a random forest on **Traitar** labels and evaluates it on **BacDive** labels, holding out one focal phylum at a time — so the model must generalize to taxa it never saw. Covers cellular traits (Oxygen_Preference, Gram_Status, Motility, Spore_Formation) for SNE (6 pretraining sizes × 5 replicates), DNABERT2 and Phylo-PCA (full size), and 8 sugar traits ("builds acid from" in BacDive vs Traitar metabolic columns) for all three embeddings. AUC is reported per trait × phylum × embedding. | `data/embedding/subset_table_*_100_*.txt`, `../../data/social_niche_embedding_100.txt`, `../../data/dnabert2_16s_embedding_reduced_100.txt`, `../../data/phylo_embed_PCA_100.txt`, `../../data/taxmap_slv_ssu_ref_nr_138.2.txt`, `data/trait_predcit.csv`, `data/bacDive.csv`, `data/agg_bac.csv` | `data/auc_res.csv` (cellular traits), `data/predict_metabolics_res.csv` (sugar traits) |
| [`traits_results.ipynb`](traits_results.ipynb) | **Main notebook** — assembles Figure 2 and Extended Data Fig. 5 from the files above (see next section). | everything listed in Section 5 | `results/*.pdf` |

## 5. The main notebook, panel by panel

### Figure 2 (four panels, `results/Figure2.pdf`)

| Panel | Title | What it shows | Data read |
|---|---|---|---|
| **A** | "facultatively − aerobic ≈ anaerobic" | Per-taxon cosine similarity between the embedding vector `v(facultatively) − v(aerobic)` and `v(anaerobic)` (Traitar labels), violin plots per group, with p values from a **phylogenetic regression** (`phylolm`, Pagel's λ) against aerobic and facultatively anaerobes. | `data/Traitar_Facultatively_Anaerobic_Anaerobic_Aerobic_all_co.csv`, `../../data/SSURefNR99_1200_slv_138_2_subset.tre` |
| **B** | R2 vs Q2 of PLS-DA | One point per metabolic trait: Traitar "Growth: Sugar" + BacDive "assimilation" / "builds acid from". Points colored red when the phylogenetic permutation pQ2 ≤ 0.05; traits on the resolution floor are ranked by z(Q2). | `data/{traitar,bacdive}/phyloperm_*.rds`, `data/traits.tsv` (Traitar accession → name map), `data/trait_predcit.csv`, `data/agg_bac.csv` |
| **C** | Trait prediction across embeddings | AUC per trait and phylum (points) with mean ± SD per embedding, for SNE (red), DNABERT-2 (blue) and Phylo-PCA (green); cellular traits block on top, sugar block shaded below. | `data/auc_res.csv` (repeat 1 at 210,000 samples), `data/predict_metabolics_res.csv` |
| **D** | Pretraining data size | AUC of the 4 cellular traits as a function of SNE pretraining corpus size (80k / 160k / 210k samples), faceted per trait, phylum as shape. | `data/auc_res.csv` (social_niche rows) |

### Extended Data Fig. 5 (`results/ExtendedDataFig5_AD.pdf`, `results/ExtendedDataFig5_EF.pdf`)

Per-trait **PLS-DA score scatter plots** (components 1 vs 2), taxa colored by
trait class, annotated with R2, Q2 and the phylogenetic p values, for:

| Panel | Source | Traits |
|---|---|---|
| A | BugBase | Gram_Status, Oxygen_Preference |
| B | Traitar | Gram_Status, Oxygen_Preference, cell_shape, Spore_Formation, Motility |
| C | BacDive | gram_stain, Oxygen_Preference, cell_shape, spore_formation, motility |
| D | Traitar | Growth: Sugar traits (ordered by Q2) |
| E | BacDive | assimilation terms |
| F | BacDive | builds_acid_from terms |

These panels reuse the same `data/{bugbase,traitar,bacdive}/phyloperm_*.rds`
files (plus `data/traits_bugbase.csv`, `data/trait_predcit.csv`,
`data/bacDive.csv`, `data/agg_bac.csv`) as Figure 2B.

## 6. Data organization

```
traits/
├── traits_annotation/            # Step 1: build the trait tables
│   ├── run_vsearch.sh            #   vsearch alignment vs BugBase reference
│   ├── Bugbase_predict.py        #   map BugBase traits onto query sequences
│   ├── Bugbase_database.txt      #   BugBase reference trait table
│   ├── traitar_predict.sh        #   Traitar phenotype prediction (SLURM)
│   ├── samples.txt               #   genome list for Traitar
│   ├── traits_predict_Bugbase.txt
│   └── traits_predict_Traitar.csv
├── data/
│   ├── 99_otus.fasta             # BugBase reference DB (Greengenes 99% OTUs)
│   ├── feces_seq_16S_silva.fasta   # query 16S sequences
│   ├── vsearch.out               # alignment output (blast6out)
│   ├── default_traits_precalculated.txt   # BugBase traits for references
│   ├── traits_precalculated.txt  # BugBase traits relabeled with query IDs
│   ├── traits_bugbase.csv        # curated BugBase table (used downstream)
│   ├── trait_predcit.csv         # curated Traitar table (used downstream)
│   ├── bacDive.csv               # BacDive table
│   ├── agg_bac.csv               # BacDive term aggregation (assimilation /
│   │                             #   builds_acid_from)
│   ├── traits.tsv                # Traitar trait accession ↔ name mapping
│   ├── embedding/
│   │   └── subset_table_{1w..16w}_100_{1..5}.txt   # subsample SNEs (Fig. 2D)
│   ├── bugbase/                  #   phyloperm_<trait>.rds  (from run_plsda.R)
│   ├── traitar/                  #   phyloperm_<trait>.rds
│   ├── bacdive/                  #   phyloperm_<trait>.rds, model_*.rds (legacy)
│   ├── logs/                     #   per-trait error logs
│   ├── auc_res.csv               # cellular-trait prediction AUCs (Fig. 2C/D)
│   ├── predict_metabolics_res.csv# sugar-trait prediction AUCs (Fig. 2C)
│   ├── fig2b_phylo_perm.csv      # flattened PLS-DA summary (Fig. 2B source)
│   ├── metabolic_res.csv         # flattened Traitar sugar summary (legacy)
│   ├── Traitar_..._all_co.csv    # analogy cosines (Fig. 2A)
│   ├── *_..._all_phy.csv         # analogy intermediates (legacy)
│   ├── pruned_tree.tre, tree_14k.nwk         # tree intermediates (legacy)
│   └── bac_KO_predicted.tsv, bac_marker_predicted_and_nsti.tsv,
│       KO_phyloglm.csv, trait_16s_output_pick.tsv   # KO/marker analyses
│                                     # from earlier iterations (not used by
│                                     # the current notebooks)
├── run_plsda.R                   # Step 2: PLS-DA + phylogenetic permutation
├── traits_predict.ipynb          # Step 3: cross-database trait prediction
├── traits_results.ipynb          # Step 4: Figure 2 + Extended Data Fig. 5
└── results/                      # Figure 2.pdf, ExtendedDataFig5_AD.pdf,
                                  #   ExtendedDataFig5_EF.pdf + intermediate
                                  #   panels (p1.pdf, p_point.pdf, merge_p*.pdf)
```

## 7. How to reproduce

### 7.1 Annotate the traits

```bash
cd traits_annotation
sbatch run_vsearch.sh                       # or run the vsearch command directly
python Bugbase_predict.py                   # BugBase traits -> traits_predict_Bugbase.txt
sbatch traitar_predict.sh                   # Traitar traits (SLURM, 28 CPUs)
```

The curated downstream tables (`data/traits_bugbase.csv`, `data/trait_predcit.csv`,
`data/bacDive.csv`) are already in place; regenerate them from the raw outputs
only if the annotations change.

### 7.2 PLS-DA with the phylogenetic null model

```bash
Rscript run_plsda.R
```

Writes one `data/<source>/phyloperm_<trait>.rds` per trait. The file currently
runs the BacDive block; the BugBase and Traitar blocks are commented out —
uncomment them to (re)generate those too. Parallelized over traits with
`doParallel`; errors per trait go to `data/logs/`.

### 7.3 Cross-database trait prediction

Run [`traits_predict.ipynb`](traits_predict.ipynb) top to bottom (Python,
`scikit-learn`). This writes `data/auc_res.csv` and
`data/predict_metabolics_res.csv`. The subsample embeddings in
`data/embedding/` must exist for the scaling part.

### 7.4 Figures

Run [`traits_results.ipynb`](traits_results.ipynb) top to bottom (R with
`tidyverse`, `ropls`-generated rds files, `phylolm`, `ape`, `phytools`,
`patchwork`). It only reads the result files above — no model fitting happens
in the notebook. Figures are written to `results/`.
