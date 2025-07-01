#!/bin/bash
#SBATCH --job-name=mafft
#SBATCH -N 1
#SBATCH -n 18
#SBATCH -p gpu
#SBATCH --mem=50
#SBATCH -o mafft.out
#SBATCH -e mafft.err

module load miniconda/4.9.2 
source activate mmseqs2_env

input="Data/pick_otu.fasta"
aligned="Data/aligned.fasta"
identity_matrix="Data/identity_matrix.txt"
distance_matrix="Data/distance_matrix.txt"
# 使用MAFFT进行自动模式多序列比对
mafft --auto --thread 18 ${input} > ${aligned}

# 使用Clustal Omega计算百分比相似性矩阵
clustalo \
  -i ${aligned} \
  --percent-id \
  --distmat-out=${identity_matrix} \
  --full \
  --force --threads 18 -t DNA

clustalo \
  -i ${aligned} \
  --distmat-out=${distance_matrix} \
  --full \
  --force --threads 18 -t DNA