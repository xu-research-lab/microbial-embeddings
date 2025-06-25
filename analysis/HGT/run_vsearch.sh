#!/bin/bash

#SBATCH --job-name=function_predict
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p cu
#SBATCH -w cu03
#SBATCH --mem=250G
#SBATCH -o out.out
#SBATCH -e out.err

module load miniconda/4.9.2
source activate qiime2-amplicon-2024.2

# barrnap --quiet --threads 4 --outseq 16S.fasta genome.fasta
queries="/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/data/feces_seq_16S_new.fasta"
database="/home/dongbiao/word_embedding_microbiome/HGT/barrnap.fna"

vsearch --usearch_global ${queries} --db ${database} --id 0.99 --blast6out ../results/vsearch_blast.out
