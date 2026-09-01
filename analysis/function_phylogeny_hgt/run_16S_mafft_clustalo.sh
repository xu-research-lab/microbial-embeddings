#!/bin/bash
#SBATCH --job-name=mafft
#SBATCH -N 1
#SBATCH -n 18
#SBATCH -p gpu
#SBATCH --mem=50
#SBATCH -o ../output/mafft.out
#SBATCH -e ../output/mafft.err

module load miniconda/4.9.2 
source activate mmseqs2_env

input="../resources/genome_mapping/data/feces_seq_16S_SLIVA.fasta"
aligned="/home/dongbiao/word_embedding_microbiome/HGT/16S/aligned.fasta"
identity_matrix="/home/dongbiao/word_embedding_microbiome/HGT/16S/identity_matrix.txt"

mafft --auto --thread 18 ${input} > ${aligned}

clustalo \
  -i ${aligned} \
  --percent-id \
  --distmat-out=${identity_matrix} \
  --full \
  --force --threads 18 -t DNA