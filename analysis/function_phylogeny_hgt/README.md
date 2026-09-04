# SNEs vs. Phylogeny, Function, and HGT

This directory contains the ecological-characterization analyses of the
[SNEs](https://github.com/xu-research-lab/microbial-embeddings) project: we
examine whether microbial **Social Niche Embeddings (SNEs)** — 100-dimensional
vectors pretrained from co-occurrence patterns in >210,000 human gut
microbiomes — recapitulate **phylogeny**, capture **functional (gene-content)
similarity**, and predict **horizontal gene transfer (HGT)** between genome
pairs.

The central document of this analysis is
**[`function_phylogeny_hgt_SNes_results.ipynb`](function_phylogeny_hgt_SNes_results.ipynb)**.
It reads the result tables produced by the scripts below and produces every
figure and statistic of the study. The only exception is
[`run_SNE_phylo_function.py`](run_SNE_phylo_function.py), which draws its
figures directly (see [Section 4](#4-scripts)). This README walks through the
notebook section by section and explains what each script does and how to
reproduce the results.

---

## 1. Scientific question

SNEs are learned purely from co-occurrence — no genomes, genes, or taxonomy are
provided during pretraining. We therefore ask how much ecological signal the
embeddings carry:

1. Does SNE similarity agree with phylogenetic distance, and does this
   agreement hold across taxonomic levels (phylum to genus)?
2. Do SNEs capture functional similarity? We test this (a) globally, by
   correlating SNE distance with distance between PICRUSt2-predicted function
   profiles, and (b) per function, by asking whether the embedding predicts the
   presence/absence of each KO while controlling for phylogeny.
3. Do SNEs predict HGT between genome pairs, and does HGT frequency decline
   with phylogenetic distance or with SNE similarity?

## 2. The data

| Resource | Path | Description |
|---|---|---|
| SNE embedding | `../../data/social_niche_embedding_100.txt` | 100-d SNE vectors for 14,093 SILVA taxa (root-level canonical input) |
| Phylogenetic embedding | `../../data/phylo_embed_PCA_100.txt` | phylogenetic embedding, PCA-reduced to 100 d |
| SILVA reference tree | `../../data/SSURefNR99_1200_slv_138_2_subset.tre` | pairwise phylogenetic distances (PhyloDM) |
| SILVA taxonomy | `../../data/taxmap_slv_ssu_ref_nr_138.2.txt` | phylum-to-genus assignments per OTU |
| PICRUSt2 function tables | `data/picrust/bac_{KO,EC,CAZY}_predicted.tsv` | predicted function profiles per OTU |
| KO annotations | `data/ko_pathway.tsv`, `data/ko_names.txt` | KO id → description; the list of KOs analyzed |
| Genome pairs | `data/genome_pairs_vsearch.txt` | pairs of genomes screened for HGT |
| 16S identity matrix | `data/identity_matrix.txt` | pairwise 16S % identity (ClustalO) |
| HGT pair table | `data/hgt.csv` | per pair: `id_1`, `id_2`, `identity` (16S % id), `cosine_co` (SNE cosine), `hgt` (candidate HGT count), `phy_dis` (phylogenetic distance) |
| Binned HGT rates | `data/hgt_plot_res.csv` | `distance`, `hgt_rate`, `group` (distance-binned summary of `hgt.csv`) |
| HGT prediction results | `data/hgt_predict_res_all.csv` | `test` (fold), `group` (feature set), `labels`, `proba` |
| Raw co-occurrence matrix | `../../script/cooccur_otuembedding/table.co` + `feature-dict.csv` | GloVe co-occurrence counts between taxon pairs ("Cooccur" feature) |

## 3. The analyses

| Analysis | Question | Method | Result |
|---|---|---|---|
| SNE vs. phylogeny | Does SNE similarity mirror phylogenetic distance? | Mean pairwise SNE cosine vs. mean phylo distance per taxonomic group; Pearson r + linear regression | `data/tax_group.csv` → notebook §1 |
| SNE vs. function (global) | Does SNE distance track functional distance? | Bray–Curtis distance on PICRUSt2 profiles vs. phylo / SNE distance; Mantel test (Pearson, 100 permutations); hexbin joint plots | `Figures/Phylo_func_*.pdf`, `Figures/SNE_func_*.pdf` |
| SNE vs. function (per KO) | Does the embedding predict presence/absence of each KO? | `phyloglm` logistic regression (`logistic_MPLE`) with phylogenetic correction; likelihood-ratio test vs. null model (df = 100), q-value FDR. In parallel, `adonis2` PERMANOVA of the SNE cosine distance on KO presence/absence | `data/phylolm/*.rds` → `data/KO_phyloglm.csv`; `data/adonis/*.rds` → `data/adonis2_res.csv` (notebook §4) |
| HGT frequency | Does HGT decline with phylo distance and with SNE similarity? | HGT rate per distance bin; LOESS + linear regression; phylo-distance bins split by SNE similarity (> 0.6 vs. < 0.6) compared with a paired t-test | `data/hgt_plot_res.csv` → notebook §2 |
| HGT prediction | Can SNE / PhyloE embeddings predict HGT between genome pairs? | 5-fold cross-validation split by genome (pairs sharing a genome never cross folds); Random Forest (500 trees) on concatenated pair embeddings; AUPRC | `data/hgt_predict_res_all.csv` → notebook §3 and §5 |

## 4. Scripts

| Script | Purpose |
|---|---|
| [`function_phylogeny_hgt_SNes_results.ipynb`](function_phylogeny_hgt_SNes_results.ipynb) | **Main results notebook (R).** Reads every table above and produces all figures of the study. |
| [`run_SNE_phylo_function.py`](run_SNE_phylo_function.py) | Global function-distance analysis (`python run_SNE_phylo_function.py KO|EC|CAZY`): samples 1,000 taxa, computes Bray–Curtis function distance vs. phylogenetic distance and vs. 1 − SNE cosine, runs Mantel tests, and **saves the figures directly** to `Figures/Phylo_func_{type}.pdf` and `Figures/SNE_func_{type}.pdf` (its results are therefore not shown in the notebook). |
| [`phylolm_subsystem_predict.R`](phylolm_subsystem_predict.R) | Per-KO phylogenetic logistic regression: fits null (`Trait ~ 1`) and full (`Trait ~` 100 SNE dimensions) `phyloglm` models and saves them as `data/phylolm/co_null_model_{KO}.rds` / `co_model_{KO}.rds` (`Rscript phylolm_subsystem_predict.R <KO>`). |
| [`run_adonis.R`](run_adonis.R) | Per-KO `adonis2` PERMANOVA of the SNE cosine distance on binarized KO presence/absence → `data/adonis/adonis2_{KO}.rds` (`Rscript run_adonis.R <KO>`). |
| [`run_blastn.sh`](run_blastn.sh) | SLURM job that detects candidate HGT between two genomes: BLASTN for every genome pair in `data/genome_pairs_vsearch.txt`, keeping hits with ≥ 99% identity and ≥ 500 bp; the count of such hits per pair becomes the `hgt` value in `data/hgt.csv`. |
| [`get_HGT_result.py`](get_HGT_result.py) | Assembles `data/hgt.csv` from the BLASTN hit counts, the 16S identity matrix, SNE cosine similarities, and phylogenetic distances. |
| [`run_HGT_predict.py`](run_HGT_predict.py) | HGT prediction: 5-fold Random Forest on four feature sets (SNE, PhyloE, Cooccur, SNE+PhyloE) → `data/hgt_predict_res_all.csv`. |
| [`run_16S_mafft_clustalo.sh`](run_16S_mafft_clustalo.sh) | SLURM job: MAFFT alignment of the 16S sequences, then ClustalO pairwise % identity matrix → `data/identity_matrix.txt`. |
| [`SNE_vs_phylogenetic_distacne.ipynb`](SNE_vs_phylogenetic_distacne.ipynb) | Computes per-taxon-group mean pairwise SNE cosine and phylo distance → `data/tax_group.csv`. |

## 5. The main notebook, section by section

[`function_phylogeny_hgt_SNes_results.ipynb`](function_phylogeny_hgt_SNes_results.ipynb)
is organized as follows:

| # | Section | What it shows |
|---|---|---|
| 1 | **The correlation between SNE similarity and phylogenetic distance** | Faceted scatter (Phylum → Genus) of mean SNE cosine vs. mean phylo distance with error bars; per-level Pearson r and regression line (`data/tax_group.csv`). |
| 2 | **The HGT frequency** | HGT rate vs. phylogenetic distance: LOESS + linear regression; phylo bins split into `Phylo SNEsim > 0.6` / `< 0.6` with a paired t-test; the SNE group plotted against a secondary axis (1 − SNESim) rescaled to the phylo-distance range (`data/hgt_plot_res.csv`). |
| 3 | **HGT predict** | AUPRC per fold (points, mean ± sd, fold lines) for SNE, PhyloE, and SNE+PhyloE, with the positive-class-proportion baseline (`data/hgt_predict_res_all.csv`, PRROC). |
| 4 | **KO vs SNEs by Adonis2 and phylolm** | Aggregates the per-KO model files: likelihood-ratio tests with q-value FDR, written to `data/KO_phyloglm.csv`; barplots of selected KOs by LR statistic (phyloglm) and F statistic (adonis2). |
| 5 | **Supplement results** | AUPRC boxplot including the Cooccur feature group, complementing §3. |

## 6. Data organization

```
Figures/                          # figures drawn directly by run_SNE_phylo_function.py
│                                 #   Phylo_func_{KO,EC,CAZY}.pdf, SNE_func_{KO,EC,CAZY}.pdf
data/
├── picrust/                      # PICRUSt2-predicted function tables (KO, EC, CAZY)
├── tax_group.csv                 # per-taxon-group mean ± sd of SNE cosine & phylo distance
├── hgt.csv                       # genome-pair table: identity, SNE cosine, HGT count, phylo distance
├── hgt_plot_res.csv              # distance-binned HGT rates per group (Phylo, SNEsim split, SNE)
├── hgt_predict_res_all.csv       # per-pair HGT prediction probabilities (all folds, all feature sets)
├── KO_phyloglm.csv               # phyloglm LRT results (p, q, LR statistic, description) per KO
├── adonis2_res.csv               # adonis2 results (R², F, p) per KO
├── ko_names.txt / ko_pathway.tsv # KO list and KO → description table
├── genome_pairs_vsearch.txt      # genome pairs screened for HGT
├── identity_matrix.txt           # pairwise 16S % identity (ClustalO)
├── genome_id_file / genome_fid.json  # genome id mappings
├── phylolm/                      # co_model_{KO}.rds, co_null_model_{KO}.rds (per KO)
├── adonis/                       # adonis2_{KO}.rds (per KO)
└── blastn_results/filtered/      # {query}.fna_vs_{subject}.fna_filtered.tsv (per pair)
```

## 7. How to reproduce

### 7.1 Environment

Python packages: `pandas`, `numpy`, `scikit-learn`, `scipy`, `scikit-bio`,
`dendropy`, `phylodm`, `seaborn`, `matplotlib`, `tqdm` (the `membed` conda
environment from the repository root, see the [main
README](../../README.md), covers these). R packages: `tidyverse`, `ggpubr`,
`pROC`, `cowplot`, `patchwork`, `ggsci`, `scales`, `ggrepel`, `PRROC`,
`qvalue`, `phylolm`, `ape`, `vegan`. The alignment steps additionally need
MAFFT, ClustalO, and BLASTN.

### 7.2 SNE vs. phylogeny

Run [`SNE_vs_phylogenetic_distacne.ipynb`](SNE_vs_phylogenetic_distacne.ipynb)
to produce `data/tax_group.csv`, then notebook §1 plots it.

### 7.3 SNE vs. function

Global comparison — the figures are written directly by the script:

```bash
python run_SNE_phylo_function.py KO
python run_SNE_phylo_function.py EC
python run_SNE_phylo_function.py CAZY
```

Per-function models — one job per KO in `data/ko_names.txt` (e.g. on a cluster):

```bash
Rscript phylolm_subsystem_predict.R K00003   # -> data/phylolm/co_model_K00003.rds, co_null_model_K00003.rds
Rscript run_adonis.R K00003                  # -> data/adonis/adonis2_K00003.rds
```

Notebook §4 then aggregates all `data/phylolm/*.rds` and `data/adonis/*.rds`
files into `data/KO_phyloglm.csv` and `data/adonis2_res.csv` and draws the
barplots.

### 7.4 HGT

```bash
# 1. Pairwise 16S identity (MAFFT + ClustalO) and per-pair BLASTN (SLURM templates)
sbatch run_16S_mafft_clustalo.sh   # -> data/identity_matrix.txt
sbatch run_blastn.sh               # -> data/blastn_results/filtered/

# 2. Assemble the pair table
python get_HGT_result.py           # -> data/hgt.csv

# 3. Predict HGT from embeddings (5-fold Random Forest)
python run_HGT_predict.py          # -> data/hgt_predict_res_all.csv
```

Notebook §2 reads `data/hgt_plot_res.csv` (a distance-binned summary of
`data/hgt.csv`), and §3/§5 read `data/hgt_predict_res_all.csv`.

### 7.5 Produce the figures

Run [`function_phylogeny_hgt_SNes_results.ipynb`](function_phylogeny_hgt_SNes_results.ipynb)
from top to bottom with an R kernel. It only reads (and in §4 aggregates) the
result files above — no model training happens in the notebook.
