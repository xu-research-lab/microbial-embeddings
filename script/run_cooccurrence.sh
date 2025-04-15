#!/bin/bash
#SBATCH --job-name=cooccurrence
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p cu
#SBATCH -w cu08
#SBATCH --mem=50G
#SBATCH -o /home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/glove_8.out
#SBATCH -e /home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/glove_8.err

module load miniconda/4.9.2
source activate jupyter_notebook

id=8
# id=$(sed "${SLURM_ARRAY_TASK_ID}q;d" matrix.txt)
data_dir="/home/dongbiao/word_embedding_microbiome/Classification_prediction/data"
glove_dir="/home/dongbiao/word_embedding_microbiome/gut_microbiome_embeddings-master/GloVe-master/build"
output_dir="/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702"
output=""
### 计算otu在样本间出现的次数
# membed dict -b $data_dir/AGP_Glove.biom -d $output_dir/feature_dict_train.txt

if [ ${id} -eq 1 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_russell_rao.bin --metric russell_rao
  membed glove-train -d feature_dict_train.txt -c AGP_russell_rao.bin \
                 -r $output/russell_rao -x 0.1 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 2 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_jaccard.bin --metric jaccard --dense --cpus 14
  membed glove-train -d feature_dict_train.txt -c AGP_jaccard.bin \
                 -r $output/jaccard -x 0.2 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 3 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_faith.bin --metric faith --dense --cpus 24
  membed glove-train -d feature_dict_train.txt -c AGP_faith.bin \
                 -r $output/faith -x 0.04 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 4 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_abundance-totalsum.bin --metric abundance --dense --normalize --cpus 24
  membed glove-train -d feature_dict_train.txt -c AGP_abundance-totalsum.bin \
                 -r $output/abundance-totalsum -x 0.0001 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 5 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_abundance-percentile.bin --metric abundance --dense --percentile --cpus 24
  membed glove-train -d feature_dict_train.txt -c AGP_abundance-percentile.bin \
                 -r $output/abundance-percentile -x 0.05 --lr 0.02 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 6 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_braycurtis_normalize.bin --metric braycurtis --dense --normalize --cpus 24
  membed glove-train -d feature_dict_train.txt -c AGP_braycurtis_normalize.bin \
                 -r $output/braycurtis-totalsum -x 0.1 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 7 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_braycurtis_percentile.bin --metric braycurtis --dense --percentile --cpus 24
  membed glove-train -d feature_dict_train.txt -c AGP_braycurtis_percentile.bin \
                 -r $output/braycurtis-percentile -x 0.1 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 8 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c $output_dir/AGP_russell_rao_weight.bin --metric glove
  membed glove-train -d $output_dir/feature_dict_train.txt -c $output_dir/AGP_russell_rao_weight.bin \
                 -r $output_dir/russell_rao_weight -x 0.1 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi
