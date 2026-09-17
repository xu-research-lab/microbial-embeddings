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

queries="data/feces_seq_16S_SLIVA.fasta" ### the OTUs sequence for social niche embedding

### getting 16S sequence from each genome
### each genome id in data/high_quality_genome.tsv
barrnap --quiet --threads 4 --outseq dir/16S_${genome_id}.fasta ${genome_id}.fasta 
cat dir/16S_${genome_id}.fasta > barrnap.fna

### vsearch global research for mapping the OTUs to genome by similarity 16S sequence
bash ../../../download_release_data.sh genome_mapping  # fetch data/barrnap.fna from the GitHub Release if missing
database="data/barrnap.fna" 
vsearch --usearch_global ${queries} --db ${database} --id 0.99 --blast6out data/vsearch_blast.out
