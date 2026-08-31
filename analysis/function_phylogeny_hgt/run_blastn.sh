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

### make reference database for each genome file
for id in `less data/genome_id_file`;do
    subject=XXX
    makeblastdb -in "$subject" -dbtype nucl -out "data/databases_blastn/$(basename "$subject" .fasta)"
done


process_pair() {  
    pair=$1  
    query=$(echo "$pair" | cut -d' ' -f1).fna  
    subject=$(echo "$pair" | cut -d' ' -f2).fna  
    
    blastn -query "XXX/${query}" \
           -db ".XXX/${subject}" \
           -outfmt "6 qseqid sseqid pident length qstart qend sstart send" \
           -perc_identity 99 -word_size 28 -num_threads 1 \
           -out "data/blastn_results/${query}_vs_${subject}.tsv"  
    
    awk -F'\t' '$4 >=500 && $3 >=99' "data/blastn_results/${query}_vs_${subject}.tsv" \
        > "data/blastn_results/filtered/${query}_vs_${subject}_filtered.tsv"  
}  

export -f process_pair

 
${path}/parallel -j 28 -a data/genome_pairs_vsearch.txt
