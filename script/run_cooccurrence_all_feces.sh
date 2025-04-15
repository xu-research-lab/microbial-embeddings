#!/bin/bash
#SBATCH --job-name=cooccurrence
#SBATCH -N 1
#SBATCH -n 14
#SBATCH -p gpu
#SBATCH -w gpu01
#SBATCH --mem=50G
#SBATCH -o glove_%a.out
#SBATCH -e glove_%a.err
#SBATCH --array=1-8

module load miniconda/4.9.2
source activate jupyter_notebook

id=$(sed "${SLURM_ARRAY_TASK_ID}q;d" matrix.txt)
biom_file="/home/dongbiao/word_embedding_microbiome/databaes/16Sdatabase/human16S_feces/human_train.biom"
glove_dir="/home/dongbiao/word_embedding_microbiome/gut_microbiome_embeddings-master/GloVe-master/build"


if [ ${id} -eq 1 ];
then
  # 计算otu在样本间出现的次数
  membed dict -b $biom_file -d feature_dict_train.txt

  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_russell_rao.bin --metric russell_rao
  membed glove-train -g $glove_dir -d feature_dict_train.txt -c ../russellrao-binary_8.bin \
                 -r russell_rao -x 0.1 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 2 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_jaccard.bin --metric jaccard --dense --cpus 14
  membed glove-train -g $glove_dir -d feature_dict_train.txt -c ../jaccard_8.bin \
                 -r jaccard -x 0.2 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 3 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_faith.bin --metric faith --dense --cpus 24
  membed glove-train -g $glove_dir -d feature_dict_train.txt -c ../faith-binary_8.bin \
                 -r faith -x 0.04 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 4 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_abundance-totalsum.bin --metric abundance --dense --normalize --cpus 24
  membed glove-train -g $glove_dir -d feature_dict_train.txt -c ../abundance-totalsum_8.bin \
                 -r abundance-totalsum -x 0.0001 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 5 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_abundance-percentile.bin --metric abundance --dense --percentile --cpus 24
  membed glove-train -g $glove_dir -d feature_dict_train.txt -c ../abundance-percentile_8.bin \
                 -r abundance-percentile -x 0.05 --lr 0.02 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 6 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_braycurtis_normalize.bin --metric braycurtis --dense --normalize --cpus 24
  membed glove-train -g $glove_dir -d feature_dict_train.txt -c ../braycurtis-totalsum_8.bin \
                 -r braycurtis-totalsum -x 0.1 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 7 ];
then
  # membed cooccur -b $data_dir/AGP_Glove.biom -c AGP_braycurtis_percentile.bin --metric braycurtis --dense --percentile --cpus 24
  membed glove-train -g $glove_dir -d feature_dict_train.txt -c ../braycurtis-percentile_8.bin \
                 -r braycurtis-percentile -x 0.1 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi

if [ ${id} -eq 8 ];
then
  membed cooccur -b $biom_file -c ../russell_rao_weight_8.bin --metric glove
  membed glove-train -g $glove_dir -d feature_dict_train.txt -c ../russell_rao_weight_8.bin \
                 -r russell_rao_weight -x 100 --lr 0.05 \
                 --embedding-size 100 --iter 100 --cpus 28
fi
