# Disease prediction

This module contains the retained classification, ablation, biomarker, and host
attribute analyses associated with Fig. 5A-I, Extended Data Fig. 7D-S, and
Extended Data Fig. 8. The classification implementation is provided by the
repository's current `membed/` package; retained wrappers and notebooks here
do not establish an adopted run.

| Subanalysis | Core files | Inputs | Outputs | Status |
| --- | --- | --- | --- | --- |
| IBD/CRC classification | `ibd_crc_classification/attention/attention_IBD_CRC.sh`, `random_forest/RF_IBD_CRC.sh`, `evaluation/TPR_FPR.R` | IBD/CRC cohort BIOMs, metadata, SNE and per-study outputs under `data/` | task figures and summaries belong under `ibd_crc_classification/results/` | Legacy source; exact run identities unknown |
| Pretraining-size ablation | `pretraining_size_ablation/IBD_CRC_datasize_attention.ipynb`, `attention/`, `pretraining/` | size-specific embeddings, cohorts and run outputs under `data/pretraining_datasize/` | `summary_ROC_results_datasize_{IBD,CRC}.csv` under `pretraining_size_ablation/results/` | Legacy source; run chain incomplete |
| Shuffled-embedding ablation | `shuffled_embedding_ablation/attention_IBD_CRC.sh` | shuffled embeddings and cohort inputs under `data/shuffle_embedding/` | compact candidate under `shuffled_embedding_ablation/results/` | Legacy source; shuffle identity unresolved |
| IBD subtype | `ibd_subtype/IBD_subtype_loo.ipynb`, `attention_UC_CD.sh` | subtype cohorts, metadata and embeddings under `data/IBD_subtype_data/` | `ibd_subtype/results/` | Legacy source; no compact output retained |
| All-disease LOO | `all_disease_loo/all_diesease_loo.ipynb`, `attention/`, `random_forest/` | disease cohorts, SNE, model checkpoints and per-study curves under `data/disease_data/` | compact summaries and figures under `all_disease_loo/results/` | Legacy source; per-study runs are not frozen |
| Biomarker analysis | `biomarker_analysis/Dimension_reduction.ipynb`, `plot_MDS.R`, `run_shap_python_{ibd,crc}.py`, `linda_all_disease.R` | cohort/model outputs and LinDA inputs under `data/` | `MDS_{IBD,CRC}_shap.csv` under `biomarker_analysis/results/` | Legacy source; LinDA and SHAP chains incomplete |
| Host attributes | `host_attribute_prediction/age/`, `host_attribute_prediction/sex/` | age/sex cohorts and metadata under `data/` | task-local `results/` directories | Legacy source; manuscript linkage unresolved |
| Shared helpers | `shared/random_forest_model/`, `shared/statistics/delong.py` | determined by callers | caller-owned results | Cross-task support; no independent result directory |

Large per-study attention/RF artifacts remain in their original relative
cohort layout under `data/` because notebooks consume checkpoints, scores,
ROC curves and images together. Compact summaries belong to the producing
subanalysis `results/`; images belong under its `results/figures/`.

This module still lacks cohort/split manifests, SNE identities for each run,
seeds/checksums, adopted aggregations, a verified LinDA generator, and a
single confirmed shuffled-embedding definition.
