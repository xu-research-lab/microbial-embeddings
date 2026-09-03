#!/bin/bash
#SBATCH --job-name=carvem
#SBATCH -N 1
#SBATCH -n 7
#SBATCH -p cu
#SBATCH --mem=50
#SBATCH -o output/carvem_%a.out
#SBATCH -e output/carvem_%a.err
#SBATCH --array=1
#SBATCH --exclude=cu01,cu02

module load miniconda/4.9.2
source activate prokka

genome_file="/home/dongbiao/word_embedding_microbiome/HGT/HGT_network/genome_id.txt"
id=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${genome_file})
# id="GCF_958301595.1_genomic"
fna_index="/home/dongbiao/word_embedding_microbiome/HGT/high_quality_genome/${id}.fna"
file_path="/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_genome"

prokka "$fna_index" --cpus 7 --outdir "${file_path}/faa/prokka_${id}" --norrna
cat "${file_path}/faa/prokka_${id}/"*.faa > "${file_path}/faa/${id}.faa"
rm -rf "${file_path}/faa/prokka_${id}"

faa_index="${file_path}/faa/${id}.faa"
mediadb="/home/dongbiao/word_embedding_microbiome/modelseed/media-main/media/WD.tsv"
model_output="/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_genome/${id}.xml"

source activate carvem_3.7
carve ${faa_index} --mediadb ${mediadb} --output ${model_output} --gapfill WD