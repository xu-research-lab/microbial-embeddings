# Disease prediction

This module contains the retained disease-classification, ablation, and
biomarker analyses for Fig. 5A-I and Extended Data Fig. 7D-S. The reusable
attention implementation is provided by the repository's `membed/` package.
Age and sex prediction are not part of the retained manuscript analysis and
are intentionally excluded.

| Analysis in manuscript order | Core code | Shared inputs | Results |
| --- | --- | --- | --- |
| IBD/CRC classification | `ibd_crc_classification/attention/attention_IBD_CRC.sh`, `random_forest/RF_IBD_CRC.sh`, `evaluation/TPR_FPR.R` | `data/IBD_CRC_model/` | `ibd_crc_classification/results/` |
| Pretraining-size ablation | `pretraining_size_ablation/IBD_CRC_datasize_attention.ipynb`, `attention/`, `pretraining/` | `data/pretraining_datasize/` | `pretraining_size_ablation/results/` |
| Shuffled-embedding ablation | `shuffled_embedding_ablation/attention_IBD_CRC.sh` | `data/shuffle_embedding/` | `shuffled_embedding_ablation/results/` |
| IBD subtype classification | `ibd_subtype/IBD_subtype_loo.ipynb`, `attention_UC_CD.sh` | `data/IBD_subtype_data/` | `ibd_subtype/results/` |
| All-disease leave-one-out | `all_disease_loo/all_disease_loo.ipynb`, `attention/`, `random_forest/` | `data/disease_data/`, `data/loo_all_*` | `all_disease_loo/results/` |
| Biomarker analysis | `biomarker_analysis/Dimension_reduction.ipynb`, `linda_all_disease.R`, `run_shap_python_{ibd,crc}.py`, `plot_MDS.R` | `data/IBD_CRC_model/`, `data/biomark/` | `biomarker_analysis/results/` |
| Shared support | `shared/random_forest_model/`, `shared/statistics/delong.py` | caller-owned | caller-owned |

`data/` is shared module input and retained run material. Each top-level
subanalysis has a `data -> ../data` link so its scripts can use that shared
layout. Nested script folders retain their own relative `data` links only
where required to run a script from that directory. Results belong to the
subanalysis that produces them; figures belong below `results/figures/`.

The retained sources do not establish a single adopted rerun. Cohort/split
manifests, per-run SNE identities, seeds/checksums, and some aggregation
records remain to be documented separately.
