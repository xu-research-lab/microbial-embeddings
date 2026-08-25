python run_attention_biom_with_SNEs.py --tasks all --dry-run

python run_attention_biom_with_SNEs.py --tasks disease --gpus 0 1 2 3 4 5 6 7 --inner-split loso \
       --run-name results_with_SNEs --linear-branch --report-combiner prob
       
python run_attention_biom_with_SNEs.py --tasks ibd_subtype --gpus 0 1 2 3 4 5 6 7 --inner-split loso \
       --run-name results_with_SNEs --linear-branch --report-combiner prob
       
python run_attention_biom_with_SNEs.py --tasks lodo --gpus 0 1 2 3 4 5 6 7 \
    --inner-split disease_loso --valid-auc macro \
    --n-estimators 32 --report-combiner prob \
    --no-linear-branch --run-name results_with_SNEs

python run_attention_biom_with_SNEs.py --tasks loso_all --gpus 0 1 2 3 4 5 6 7 \
    --inner-split same_disease --set loss=LogitAdjusted --logit-adjust-tau 1 \
    --n-estimators 1 --report-combiner prob \
    --run-name results_with_SNEs --cv-folds 5 --no-linear-branch

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
