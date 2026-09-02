# Disease Classification with Leave-One-Out Validation

This directory contains the disease-classification benchmark of the
[SNEs](https://github.com/xu-research-lab/microbial-embeddings) project: we evaluate
whether microbial **Social Niche Embeddings (SNEs)** — 100-dimensional vectors
pretrained from co-occurrence patterns in >210,000 human gut microbiomes — improve
deep-learning-based disease classification compared to classical baselines
(Random Forest, SVM) and to alternative microbial representations (phylogenetic
embeddings, DNABERT2).

The central document of this analysis is
**[`disease_classification.ipynb`](disease_classification.ipynb)**. It reads the
results produced by the scripts below and produces every figure and statistic of
the study. This README walks through the notebook section by section and explains
what each script does and how to reproduce the results.

---

## 1. Scientific question

Gut microbiome studies are notoriously hard to generalize across cohorts
(batch effects, geography, sequencing platforms). We ask:

1. Does an attention model operating on SNE-embedded taxa classify diseases
   better than Random Forest / SVM on raw abundances, when every evaluation
   holds out an **entire study** that the model never saw during training?
2. Does this hold when a **whole disease** is held out?
3. Are the embeddings actually responsible for the performance (shuffle controls,
   alternative embeddings, pretraining-data-size scaling)?
4. Which taxa drive the model's decisions (attention, integrated gradients,
   SHAP), and do they match established biomarkers (LinDA differential
   abundance)?

## 2. The data

| Resource | Path | Description |
|---|---|---|
| Pretraining BIOM table | `../../data/gut_pretraining.biom` | ~210,090 samples x 14,093 SILVA taxa, used to pretrain the SNEs |
| SNE embedding (main) | `../../data/social_niche_embedding_removing_disease_samples_100.txt` | 100-d SNE vectors (disease samples removed during pretraining) |
| Shuffled SNE (control) | `../../data/social_niche_embedding_removing_disease_samples_100_shuffled.txt` | column-wise shuffled SNE matrix |
| Phylogenetic embedding | `../../data/phylo_embed_PCA_100.txt` | phylogenetic embedding, PCA-reduced to 100 d |
| DNABERT2 embedding | `../../data/dnabert2_16s_embedding_reduced_100.txt` | DNABERT2 16S embedding, reduced to 100 d |
| Sample metadata | `../../data/metadata_disease_classification.tsv` | study, disease (`disease_name_ab`), control/case (`group`), etc. |

Thirteen diseases are analyzed: **AS, ASD, BD, CAD, CRC, GD, IBD, IBS, MS, OB, PD, SZ, T2DM**
(57 disease–study cohorts in total, see
[`run_leave_one_study_out_each_diease_list.tsv`](run_leave_one_study_out_each_diease_list.tsv)).
A cohort is kept only if it has >= 25 samples per group and contains both classes
(notebook cells 6–7).

## 3. The models

| Method | Meaning |
|---|---|
| `RF` | Random Forest on normalized raw abundances |
| `SVM` | SVM on normalized raw abundances |
| `RF_SNEs` | Random Forest on SNE features |
| `RF_PCA` | Random Forest on phylogenetic-PCA features |
| `Att_SNEs` | Attention model (`membed` `OtuAttentionEncoder`) with SNE embeddings |
| `Att_SNEs_shuffle` | Attention model with the shuffled SNE matrix (negative control) |
| `Att_both_shuffle` | Attention model with shuffled SNE **and** globally shuffled abundance table |
| `Att_PhyloE` | Attention model with phylogenetic embeddings |
| `Att_DNABERT2` | Attention model with DNABERT2 embeddings |

### Validation schemes

| Task (runner) | TSV defining folds | Held out | Result location |
|---|---|---|---|
| `disease` | `run_leave_one_study_out_each_diease_list.tsv` | one study of one disease (`{disease}_{study}`) | `Data/disease_data/` |
| `ibd_subtype` | `run_leave_one_study_out_IBD_subtype_list.tsv` | one study of one IBD subtype (UC, CD, CCD, ICD, feces/biopsy splits) | `Data/IBD_subtype_data/` |
| `lodo` | `run_leave_one_disease_out.tsv` | one whole disease | `Data/loo_all_diseases/` |
| `loso_all` | `run_leave_one_study_out.tsv` | one study in the pooled all-disease table | `Data/loo_all_studies/` |
| `shuffled_table` | `run_shuffled_table.tsv` | like `disease`, but on globally shuffled BIOM tables | `Data/shuffle_table_IBD_CRC/` |

## 4. Scripts

| Script | Purpose |
|---|---|
| [`run_attention_biom_with_SNEs.py`](run_attention_biom_with_SNEs.py) | **Main trainer.** Runs the Attention model for all tasks. |
| [`run_rf.py`](run_rf.py) | Random Forest baseline: rank-normalize. |
| [`run_svm.py`](run_svm.py) | SVM baseline with the same preprocessing. |
| [`explain_attention.py`](explain_attention.py) | **Model interpretation.** Per-taxon attributions for every fold: attention readout (`attn`, reported as enrichment) and integrated gradients (`ig`). |
| [`get_subdatasets.py`](get_subdatasets.py) | Subsamples the pretraining table (1k–160k samples, 5 replicates) for the scaling experiment. |
| [`SNEs_shuffle.py`](SNEs_shuffle.py) | Shuffles each embedding column independently (row labels kept) -> shuffled-SNE control. |
| [`biom_table_shuffle.py`](biom_table_shuffle.py) | Globally permutes all values of a BIOM table -> shuffled-abundance control. |
| [`run_jobs.sh`](run_jobs.sh) | **Master script** with the exact commands used for every experiment in the study (see reproduction below). |

### Run-list TSVs

The `run_*.tsv` files are the fold definitions; each row is one fold
(disease/study pairs). The runner reads them via the `TASKS` table in
`run_attention_biom_with_SNEs.py:106`.

## 5. The main notebook, section by section

[`disease_classification.ipynb`](disease_classification.ipynb) is organized as
follows:

| # | Section | What it shows |
|---|---|---|
| 1 | **AUC Performance of RF and Attention Models in IBD and CRC** | Per-study AUC boxplots comparing RF, SVM, RF_SNEs, RF_PCA, Att_SNEs, Att_SNEs_shuffle, Att_PhyloE, Att_DNABERT2. |
| 2 | **ROC plots: SNEs vs RF (IBD, CRC)** | One ROC curve per study (solid = SNEs, dashed = RF) plus a mean curve, with an AUC/count table. |
| 3 | **Pretraining datasize** | AUC as a function of pretraining corpus size (1k–160k samples) for IBD and CRC. |
| 4 | **Shuffle SNEs and biom table** | AUC degradation across `SNE -> Shuffled SNE -> Shuffled Both` (negative controls). |
| 5 | **Leave one dataset out** | Per-study AUC of RF vs SNEs across all studies, paired t-test. |
| 6 | **Leave one disease out** | Per-study "macro" AUC per left-out disease, paired tests, slope plot. |
| 7 | **Model explain (grad)** | MDS of biomarker taxa in SNE vs phylo space, silhouette scores with 999 label permutations. |
| S1 | **IBD and CRC model attention weight** | Mean attention enrichment of LinDA-defined control/disease markers vs non-markers (Kruskal-Wallis + Mann-Whitney U). |
| S2 | **Attention heatmap** | Samples x pooling-dimension heatmaps of the pooled encoder representation, hierarchically clustered, checked for fold-driven (batch) splits. |
| S3 | **MDS** | MDS of biomarker groups per disease. |
| S4 | **IBD subtype** | 12-panel ROC grid (UC/CD/CCD/ICD x feces/biopsy/feces_to_biopsy), SNEs vs RF. |
| S5 | **ROC: SNEs vs RF, other diseases** | ROC grids for the remaining 11 diseases (AS, ASD, BD, CAD, GD, IBS, MS, OB, PD, SZ, T2DM). |

## 6. Companion notebooks

| Notebook | Purpose |
|---|---|
| [`Biomarker_analysis.ipynb`](Biomarker_analysis.ipynb) | In-depth biomarker analysis of the `explain_attention.py` outputs: marker selection, MDS + PERMANOVA, silhouette permutation tests, attention vs LinDA markers, pooled heatmaps. |
| [`disease_classification_embedding.ipynb`](disease_classification_embedding.ipynb) | Early exploration: weighted-SNE features and RF results for CRC/IBD. |
| [`PCA_reduction.ipynb`](PCA_reduction.ipynb) | PCA reduction of the phylogenetic embedding to 100 d. |

## 7. Data organization

```
Data/
├── disease_data/                 # 'disease' task (per-disease LOSO)
│   ├── {disease}/                #   e.g. CRC/
│   │   ├── {study}/              #     PRJEB6070/
│   │   │   ├── train_loo.biom / test_loo.biom
│   │   │   ├── results/          # attention_scores.csv, roc_curve.csv,
│   │   │   │                     # auc_loo.png, loss_loo.png, attention_loo.pt
│   │   │   └── RF/               # RF_Scores.csv, RF_ROC.csv
│   │   └── metadata.tsv
│   ├── results_with_SNEs/        # one dir per fold, members/*/pred_test.csv
│   ├── results_with_shuffled_SNEs/
│   ├── results_with_phylo_embed_PCA/
│   ├── results_with_dnabert2/
│   ├── _results_with_rf/         # rf/pred_test.csv per fold
│   ├── _results_with_svm/        # svm/pred_test.csv per fold
│   ├── results_with_rf_embed/    # RF on SNE features
│   └── results_with_rf_PCA/      # RF on phylo-PCA features
├── loo_all_studies/              # 'loso_all' task (results_with_SNEs, _results_with_rf)
├── loo_all_diseases/             # 'lodo' task (pooled all-disease table + results)
├── IBD_subtype_data/             # 'ibd_subtype' task
├── IBD_CRC_model/                # LinDA differential-abundance results per disease
│                                 #   (linda_res_study.csv, model-independent markers)
├── shuffle_table_IBD_CRC/        # 'shuffled_table' task
├── pretraining_datasize/         # scaling experiment
│   ├── subset/                   #   subsampled tables (get_subdatasets.py)
│   ├── embedding/                #   SNEs re-trained per size (run_cooccur_SNEs.sh)
│   ├── trainning_data/           #   LOSO folds per size
│   └── result/                   #   summary_ROC_results_datasize_{IBD,CRC}.csv
└── biomark/                      # explain_attention.py outputs
    ├── attn_summary_{task}.csv          # fold -> disease, attention readout
    ├── shap_grad_{task}_{fold}.csv      # per-sample per-taxon attributions
    ├── shap_summary_{task}.csv
    ├── pool_{task}_{fold}.csv           # pooled encoder representation per sample
    ├── MDS_*.csv                        # MDS coordinates of marker taxa
    └── attn/shap_consistency_*.csv      # cross-cohort consistency tables
```

Per-fold prediction files (`pred_test.csv`) contain `sample_id`, `true_label`,
`prob`; ensemble members live in `members/*/pred_test.csv` and are averaged by
the notebook's `load_prediction()`.

## 8. How to reproduce

### 8.1 Environment

From the repository root (`microbial-embeddings/`):

```bash
conda env create --name membed --file requirements_dev.yml
conda activate membed
pip install -e .
cd analysis/Disease_classification_loo
```

Training uses GPUs (the commands below assume 8; adjust `--gpus`).

### 8.2 Run the classification experiments

[`run_jobs.sh`](run_jobs.sh) contains the exact commands used. In short:

```bash
# Dry-run: list all jobs without running anything
python run_attention_biom_with_SNEs.py --tasks all --dry-run

# Main per-disease LOSO run with SNE embeddings (Attention model)
python run_attention_biom_with_SNEs.py --tasks disease --gpus 0 1 2 3 4 5 6 7 \
    --inner-split loso --run-name results_with_SNEs \
    --linear-branch --report-combiner prob

# IBD subtypes
python run_attention_biom_with_SNEs.py --tasks ibd_subtype --gpus 0 1 2 3 4 5 6 7 \
    --inner-split loso --run-name results_with_SNEs \
    --linear-branch --report-combiner prob

# Leave-one-disease-out and leave-one-study-out (pooled table)
python run_attention_biom_with_SNEs.py --tasks lodo --gpus 0 1 2 3 4 5 6 7 \
    --inner-split per_disease --set loss=GroupBalanced --group-balance-beta 0.5 \
    --logit-adjust-tau 1 --no-linear-branch \
    --report-combiner prob --run-name results_with_SNEs

python run_attention_biom_with_SNEs.py --tasks loso_all --gpus 0 1 2 3 4 5 6 7 \
    --inner-split per_disease --set loss=LogitAdjusted --logit-adjust-tau 1 \
    --report-combiner prob --run-name results_with_SNEs --no-linear-branch

# Alternative embeddings and controls (same task, different --glove-embedding)
python run_attention_biom_with_SNEs.py --tasks disease --gpus 0 1 2 3 4 5 6 7 \
    --run-name results_with_dnabert2 --linear-branch --report-combiner prob \
    --inner-split loso \
    --glove-embedding ../../data/dnabert2_16s_embedding_reduced_100.txt

python run_attention_biom_with_SNEs.py --tasks disease --gpus 0 1 2 3 4 5 6 7 \
    --run-name results_with_phylo_embed_PCA --linear-branch --report-combiner prob \
    --inner-split loso --glove-embedding ../../data/phylo_embed_PCA_100.txt

python run_attention_biom_with_SNEs.py --tasks disease --gpus 0 1 2 3 4 5 6 7 \
    --run-name results_with_shuffled_SNEs --linear-branch --report-combiner prob \
    --inner-split loso \
    --glove-embedding ../../data/social_niche_embedding_removing_disease_samples_100_shuffled.txt

# Shuffled-abundance control (create the tables first)
python biom_table_shuffle.py
python run_attention_biom_with_SNEs.py --tasks shuffled_table --gpus 0 1 2 3 4 5 6 7 \
    --run-name results_with_shuffled_both --linear-branch --report-combiner prob \
    --inner-split loso \
    --glove-embedding ../../data/social_niche_embedding_removing_disease_samples_100_shuffled.txt

# Baselines
python run_rf.py  --tasks all
python run_svm.py --tasks all
```

### 8.3 Model interpretation

```bash
# IBD & CRC attention / gradients / SHAP
python explain_attention.py --task disease --run-name results_with_SNEs_ckpt \
    --run-tsv run_leave_one_study_out_explain.tsv \
    --diseases CRC IBD --gpus 0 1 2 3 4 5 6 7

# Leave-one-disease-out attributions
python explain_attention.py --task lodo --run-name results_with_SNEs_ckpt \
    --gpus 0 1 2 3 4 5 6 7 --perm-per-disease 0
```

`explain_attention.py` writes three readouts per fold into `Data/biomark/`:
- `attn` — mean attention distribution ("which microbes does the model attend
  to"), reported as enrichment where 1.0 = an average taxon; also records the
  attention-vs-abundance Spearman correlation per fold;
- `ig` — integrated gradients;

### 8.4 Produce the figures

Run [`disease_classification.ipynb`](disease_classification.ipynb) from top to
bottom (jupyter/lab, `membed` environment). It only reads the result files above
— no training happens in the notebook. Figures are written to `Figures/` and
`Data/biomark/`.
