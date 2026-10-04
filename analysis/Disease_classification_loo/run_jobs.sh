# Commands behind the results under Data/ that the notebooks read.
# Every flag below matches what the result.json / split.json files on disk
# record; where a run was later renamed, the original name is noted.

### >>> 2026-10-03: PRJNA578223 (ASD) DROPPED -- its ASD and control groups were sequenced
### with different primers (338F/806R + barcodes vs 341F/805R). ASD = 2 studies; disease task
### 63 folds, loso_all 61 folds, lodo 13; 10,518 samples. Every run below was redone on this
### data on 2026-10-03 07:46-10:03 (gpu01 ~/tmp_claude/drop578223_runs.sh: RF/SVM --overwrite;
### disease SNEs/dnabert2/phylo with --resume = only the ASD folds; lodo and loso_all in full).
### Backup: ~/tmp_claude/backup_drop578223_20261003/  <<<
###
### >>> 2026-10-02 data update: (a)-(d) below were RUN on 2026-10-02 16:22-19:01 <<<
### (gpu01 chain ~/tmp_claude/snes_rerun_20261002.sh, --workers-per-gpu 2; SVM on cu01).
### (e) was not run. The mkdir/mv/svm lines are commented out: rerunning them would move
### the new results away.
# Data now: disease task 64 folds (GD = 4 studies, PRJNA1250469 dropped; IBS = 7 public
# cohorts, DADA2 reads truncated to 150 bp, PRJEB44533 121 bp). lodo (13 folds) and
# loso_all (62 folds: +7 IBS, -PRJNA1250469, -PRJNA268708) were rebuilt from the disease
# task as union profiles, test features NOT aligned to the training axis;
# run_leave_one_study_out.tsv has 62 rows. Old data/results:
# ~/tmp_claude/backup_unaligned_switch_20261002/
# (a) lodo, all 13 folds:      section 2 command, unchanged
# (b) loso_all, all 62 folds:  section 2 command, unchanged
# (c) section 3, results_with_dnabert2 and results_with_phylo_embed_PCA date from 08-25:
#     they lack the 7 IBS folds and GD_PRJNA450230/GD_PRJNA799831, and their GD, SZ and
#     T2DM folds are stale (data changed 09-30). Move them aside, then run both section-3
#     commands (all 64 folds):
# mkdir -p ~/tmp_claude/backup_disease_embed_20261002
# mv Data/disease_data/results_with_dnabert2 Data/disease_data/results_with_phylo_embed_PCA \
#    Data/disease_data/_results_with_svm ~/tmp_claude/backup_disease_embed_20261002/
# (d) section 5 SVM (CPU) has the same gaps; after the move above:
# python run_svm.py --tasks disease --run-name results_with_svm
# (e) only if the lodo attention figures are used: section 6 lodo checkpoint run, then
#     explain_attention.py --task lodo (the lodo data changed).
# Up to date, no rerun: section 1 (disease SNEs and RF), ibd_subtype, the CRC/IBD-only
# runs in sections 3 and 6, section 4, RF for disease/lodo/loso_all.

python run_attention_biom_with_SNEs.py --tasks all --dry-run

### 0. controls and data prep
python SNEs_shuffle.py --seed 5      # -> ../../data/..._100_shuffled.txt
python biom_table_shuffle.py         # -> Data/shuffle_table_IBD_CRC/
python get_subdatasets.py            # -> Data/pretraining_datasize/trainning_data/
sbatch run_cooccur_SNEs.sh           # -> Data/pretraining_datasize/embedding/

### 1. single-disease LOSO (Data/disease_data/results_with_SNEs, 63 folds since 2026-10-03;
###    GD rerun 2026-10-01, IBS 2026-10-02, ASD 2026-10-03)
python run_attention_biom_with_SNEs.py --tasks disease --gpus 0 1 2 3 4 5 6 7 --inner-split loso \
       --run-name results_with_SNEs --report-combiner prob

# NOTE: Extended Data Fig. 7 (IBD subtypes) still reads the older single-model
# results in Data/IBD_subtype_data/<disease>/<study>/results/; this run has
# not been made yet (Data/IBD_subtype_data/results_with_SNEs does not exist).
python run_attention_biom_with_SNEs.py --tasks ibd_subtype --gpus 0 1 2 3 4 5 6 7 --inner-split loso \
       --run-name results_with_SNEs --report-combiner prob

### 2. all diseases pooled
# leave-one-disease-out (Data/loo_all_diseases/results_with_SNEs, 13 x 12 members;
# originally run as results_with_SNEs_test_remove_lowsample). Data since 2026-10-02:
# union profiles, test features not aligned to the training axis.
python run_attention_biom_with_SNEs.py --tasks lodo --gpus 0 1 2 3 4 5 6 7 \
    --patience 2 \
    --inner-split per_disease --report-combiner prob \
    --set loss=GroupBalanced+LogitAdjusted --group-balance-beta 0.5 --logit-adjust-tau 1 \
    --run-name results_with_SNEs --no-linear-branch 

# leave-one-study-out (Data/loo_all_studies/results_with_SNEs; 701 members on the old
# 56 folds, originally run as results_with_SNEs_per_disease). Data since 2026-10-03:
# 61 folds, built like lodo (not aligned).
# --patience 2 (default 15): POOLED TASKS ONLY (lodo above and loso_all below). Adopted 2026-10-04.
# A pooled fold trains on ~10,000 samples drawn from twelve diseases, and its members overfit: the ones
# that stopped latest and reached the highest training AUC transferred worst. Cutting patience raised
# loso_all 0.677 -> 0.683 (34 folds better / 21 worse, Wilcoxon P=0.039) and lodo 0.635 -> 0.638, while
# the loso member mean AUC ROSE 0.664 -> 0.673 and the longest run fell from 65 epochs to 14.
# The single-disease task keeps the default patience 15 on purpose: its folds train on 50-250 samples,
# where patience 2 underfits -- it cost 0.707 -> 0.696 overall and far more on the folds with one training
# cohort (AS -0.074, SZ -0.066, CAD -0.035), and it erased the SNEs advantage over phylo-PCA (P 0.055 ->
# 0.474) while leaving the dnabert2/phylo baselines unchanged. Measurements: ~/tmp_claude/pat2_compare.tsv
# and the backups backup_pat15_allsnes_20261004/ (patience-15 disease results, restored) and
# backup_disease_pat2_20261004/ (the patience-2 disease run, parked).
# NOT touched: the shuffled-embedding control (section 3, last run 2026-09-17) and ibd_subtype, both of
# which predate the current data.
python run_attention_biom_with_SNEs.py --tasks loso_all --gpus 0 1 2 3 4 5 6 7 \
    --patience 2 \
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
python run_rf.py  --tasks disease lodo loso_all --run-name results_with_rf   # current: lodo/loso_all rerun 2026-10-02 (--overwrite)
python run_svm.py --tasks disease --run-name results_with_svm
python run_rf_CRC_continent.py
python run_rf_CRC_sample_efficiency.py

### 6. model explanation (needs checkpoints, so retrain with --keep-ckpt)
python run_attention_biom_with_SNEs.py --tasks disease \
    --gpus 0 1 2 3 4 5 6 7 --inner-split loso --linear-branch \
    --report-combiner prob \
    --run-tsv run_leave_one_study_out_explain.tsv \
    --run-name results_with_SNEs_ckpt --keep-ckpt

# lodo, RERUN 2026-10-04 21:22-21:35 (gpu01, ~/tmp_claude/lodo_explain.sh, log
# lodo_explain.log). The files it replaced dated from 09-01 and covered 13 folds
# including ASD, on data rebuilt 10-03 22:21 -- Figure 5 panel I was drawn from
# models that no longer exist. Old biomark files: ~/tmp_claude/backup_biomark_lodo_20261004/
# (its asd_leftovers/ holds the ASD ones).
#
# The flags below are section 2's lodo command plus --keep-ckpt, NOT the
# 32-member --inner-split disease_loso command used until 10-04: that one
# explained a different ensemble from the one reported, which
# explain_attention.py's own docstring flags and gives this command for.
# Verified: all 12 folds reproduce lodo.results_with_SNEs.csv to 0.000000
# (mean 0.6380) and explain_attention reports pred_diff ~1e-16 per fold.
# The 09-01 32-member checkpoints are still in results_with_SNEs_ckpt/ and are
# now stale for lodo; Data/disease_data/results_with_SNEs_ckpt/ (CRC+IBD) is not.
python run_attention_biom_with_SNEs.py --tasks lodo --gpus 0 1 2 3 4 5 6 7 \
    --patience 2 \
    --inner-split per_disease --report-combiner prob \
    --set loss=GroupBalanced+LogitAdjusted --group-balance-beta 0.5 --logit-adjust-tau 1 \
    --no-linear-branch \
    --run-name results_with_SNEs_ckpt_per_disease --keep-ckpt

## IBD, CRC
python explain_attention.py --task disease --run-name results_with_SNEs_ckpt \
    --run-tsv run_leave_one_study_out_explain.tsv \
    --diseases CRC IBD --gpus 0 1 2 3 4 5 6 7 --dump-pooled

## permutation check on one fold (Data/biomark_perm)
python explain_attention.py --task disease --run-name results_with_SNEs_ckpt \
    --run-tsv run_leave_one_study_out_explain.tsv --diseases CRC \
    --gpus 0 --perm-per-disease 1 --perm-max-folds 1 --out-dir Data/biomark_perm

## lodo (run-name follows the retrain above)
python explain_attention.py --task lodo --run-name results_with_SNEs_ckpt_per_disease \
    --gpus 0 1 2 3 4 5 6 7 --perm-per-disease 0
