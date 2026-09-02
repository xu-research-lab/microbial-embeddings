python run_attention_biom_with_SNEs.py --tasks all --dry-run

python run_attention_biom_with_SNEs.py --tasks disease --gpus 0 1 2 3 4 5 6 7 --inner-split loso \
       --run-name results_with_SNEs --linear-branch --report-combiner prob 
       
python run_attention_biom_with_SNEs.py --tasks ibd_subtype --gpus 0 1 2 3 4 5 6 7 --inner-split loso \
       --run-name results_with_SNEs --linear-branch --report-combiner prob 
       
python run_attention_biom_with_SNEs.py --tasks lodo --gpus 0 1 2 3 4 5 6 7 \
    --inner-split per_disease --report-combiner prob \
    --set loss=GroupBalanced --group-balance-beta 0.5 --logit-adjust-tau 1 \
    --no-linear-branch --run-name results_with_SNEs

python run_attention_biom_with_SNEs.py --tasks loso_all --gpus 0 1 2 3 4 5 6 7 \
    --inner-split per_disease --set loss=LogitAdjusted --logit-adjust-tau 1 \
    --report-combiner prob --run-name results_with_SNEs --no-linear-branch

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

python run_attention_biom_with_SNEs.py --tasks disease \
       --gpus 0 1 2 3 4 5 6 7 \
       --run-name results_with_shuffled_SNEs \
       --linear-branch --report-combiner prob --inner-split loso \
       --glove-embedding ../../data/social_niche_embedding_removing_disease_samples_100_shuffled.txt

python run_attention_biom_with_SNEs.py --tasks shuffled_table \
       --gpus 0 1 2 3 4 5 6 7 \
       --run-name results_with_shuffled_both \
       --linear-branch --report-combiner prob --inner-split loso \
       --glove-embedding ../../data/social_niche_embedding_removing_disease_samples_100_shuffled.txt

#shape model explain
## 1. IBD, CRC
python explain_attention.py --task disease --run-name results_with_SNEs_ckpt \
    --run-tsv run_leave_one_study_out_explain.tsv \
    --diseases CRC IBD --gpus 0 1 2 3 4 5 6 7

### 2. lodo
python explain_attention.py --task lodo --run-name results_with_SNEs_ckpt \
    --gpus 0 1 2 3 4 5 6 7 --perm-per-disease 0
