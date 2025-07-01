#!/bin/bash
#SBATCH --job-name=blastn
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p cu
#SBATCH --mem=250
#SBATCH -o blastn_%a.out
#SBATCH -e blastn_%a.err
#SBATCH --exclude=cu01

module load miniconda/4.9.2 
source activate qiime2-amplicon-2024.2

path=/home/dongbiao/software/parallel-20250222/bin

process_pair() {  
    pair=$1  
    query=$(echo "$pair" | cut -d' ' -f1).fna  
    subject=$(echo "$pair" | cut -d' ' -f2).fna  
    
    blastn -query "../high_quality_genome/${query}" \
           -db "../HGT_network/databases/${subject}" \
           -outfmt "6 qseqid sseqid pident length qstart qend sstart send" \
           -perc_identity 99 -word_size 28 -num_threads 1 \
           -out "../HGT_network/results/${query}_vs_${subject}.tsv"  
    
    awk -F'\t' '$4 >=500 && $3 >=99' "../HGT_network/results/${query}_vs_${subject}.tsv" \
        > "../HGT_network/results/filtered/${query}_vs_${subject}_filtered.tsv"  
}  

export -f process_pair

 
${path}/parallel -j 28 -a Data/genome_pairs_vsearch.txt
