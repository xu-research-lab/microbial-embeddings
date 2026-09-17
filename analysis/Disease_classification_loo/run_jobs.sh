# Commands behind the results under Data/ that the notebooks read.
# Every flag below matches what the result.json / split.json files on disk
# record; where a run was later renamed, the original name is noted.

python run_attention_biom_with_SNEs.py --tasks all --dry-run

### 0. controls and data prep
python SNEs_shuffle.py --seed 5      # -> ../../data/..._100_shuffled.txt
python biom_table_shuffle.py         # -> Data/shuffle_table_IBD_CRC/
python get_subdatasets.py            # -> Data/pretraining_datasize/trainning_data/
sbatch run_cooccur_SNEs.sh           # -> Data/pretraining_datasize/embedding/

### 1. single-disease LOSO (Data/disease_data/results_with_SNEs, 314 members)
python run_attention_biom_with_SNEs.py --tasks disease --gpus 0 1 2 3 4 5 6 7 --inner-split loso \
       --run-name results_with_SNEs --linear-branch --report-combiner prob

# NOTE: Extended Data Fig. 7 (IBD subtypes) still reads the older single-model
# results in Data/IBD_subtype_data/<disease>/<study>/results/; this run has
# not been made yet (Data/IBD_subtype_data/results_with_SNEs does not exist).
python run_attention_biom_with_SNEs.py --tasks ibd_subtype --gpus 0 1 2 3 4 5 6 7 --inner-split loso \
       --run-name results_with_SNEs --linear-branch --report-combiner prob

### 2. all diseases pooled
# leave-one-disease-out (Data/loo_all_diseases/results_with_SNEs, 13 x 12 members;
# originally run as results_with_SNEs_test_remove_lowsample)
python run_attention_biom_with_SNEs.py --tasks lodo --gpus 0 1 2 3 4 5 6 7 \
    --inner-split per_disease --report-combiner prob \
    --set loss=GroupBalanced+LogitAdjusted --group-balance-beta 0.5 --logit-adjust-tau 1 \
    --no-linear-branch --run-name results_with_SNEs

# leave-one-study-out (Data/loo_all_studies/results_with_SNEs, 701 members;
# originally run as results_with_SNEs_per_disease)
python run_attention_biom_with_SNEs.py --tasks loso_all --gpus 0 1 2 3 4 5 6 7 \
    --inner-split per_disease \
    --set loss=GroupBalanced+LogitAdjusted --group-balance-beta 0.5 --logit-adjust-tau 1 \
    --n-estimators 1 --report-combiner prob \
    --no-linear-branch --run-name results_with_SNEs

### 3. other embeddings and shuffled controls
python run_attention_biom_with_SNEs.py --tasks disease \
       --gpus 0 1 2 3 4 5 6 7 \
       --run-name results_with_dnabert2 \
       --linear-branch --report-combiner prob --inner-split loso \
       --glove-embedding ../../data/dnabert2_16s_embedding_reduced_100.txt

python run_attention_biom_with_SNEs.py --tasks disease \
       --gpus 0 1 2 3 4 5 6 7 \
       --run-name results_with_phylo_embed_PCA \
       --linear-branch --report-combiner prob --inner-split loso \
       --glove-embedding ../../data/phylo_embed_PCA_100.txt

# CRC and IBD folds only (98 members; originally results_with_SNE_shuffle_seed5)
python run_attention_biom_with_SNEs.py --tasks disease \
       --run-tsv run_shuffled_table.tsv \
       --gpus 0 1 2 3 4 5 6 7 \
       --run-name results_with_shuffled_SNEs \
       --linear-branch --report-combiner prob --inner-split loso \
       --glove-embedding ../../data/social_niche_embedding_removing_disease_samples_100_shuffled.txt

python run_attention_biom_with_SNEs.py --tasks shuffled_table \
       --gpus 0 1 2 3 4 5 6 7 \
       --run-name results_with_shuffled_both \
       --linear-branch --report-combiner prob --inner-split loso \
       --glove-embedding ../../data/social_niche_embedding_removing_disease_samples_100_shuffled.txt

### 4. CRC geography and sample efficiency
python run_attention_biom_CRC_continent.py --gpus 0 1 2 3 4 5 6 7 \
       --run-name crc_geo --linear-branch --report-combiner prob
python run_attention_biom_CRC_sample_efficiency.py --gpus 0 1 2 3 4 5 6 7 \
       --run-name crc_se --linear-branch --report-combiner prob

### 5. baselines (run_rf.py writes to Data/<task>/_<run-name>/)
python run_rf.py  --tasks disease lodo loso_all --run-name results_with_rf
python run_svm.py --tasks disease --run-name results_with_svm
python run_rf_CRC_continent.py
python run_rf_CRC_sample_efficiency.py

### 6. model explanation (needs checkpoints, so retrain with --keep-ckpt)
python run_attention_biom_with_SNEs.py --tasks disease \
    --gpus 0 1 2 3 4 5 6 7 --inner-split loso --linear-branch \
    --report-combiner prob \
    --run-tsv run_leave_one_study_out_explain.tsv \
    --run-name results_with_SNEs_ckpt --keep-ckpt

# The existing lodo checkpoints are this 32-member disease_loso ensemble, NOT
# the per_disease ensemble scored in section 2 -- see explain_attention.py.
python run_attention_biom_with_SNEs.py --tasks lodo \
    --gpus 0 1 2 3 4 5 6 7 --inner-split disease_loso --valid-auc macro \
    --n-estimators 32 --report-combiner prob --no-linear-branch \
    --run-name results_with_SNEs_ckpt --keep-ckpt

## IBD, CRC
python explain_attention.py --task disease --run-name results_with_SNEs_ckpt \
    --run-tsv run_leave_one_study_out_explain.tsv \
    --diseases CRC IBD --gpus 0 1 2 3 4 5 6 7 --dump-pooled

## permutation check on one fold (Data/biomark_perm)
python explain_attention.py --task disease --run-name results_with_SNEs_ckpt \
    --run-tsv run_leave_one_study_out_explain.tsv --diseases CRC \
    --gpus 0 --perm-per-disease 1 --perm-max-folds 1 --out-dir Data/biomark_perm

## lodo
python explain_attention.py --task lodo --run-name results_with_SNEs_ckpt \
    --gpus 0 1 2 3 4 5 6 7 --perm-per-disease 0
