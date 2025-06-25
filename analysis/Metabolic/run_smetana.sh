#!/bin/bash
#SBATCH --job-name=smetana
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p cu
#SBATCH --mem=100
#SBATCH -o ../output/smetana_%a.out
#SBATCH -e ../output/smetana_%a.err
#SBATCH --array=19
#SBATCH --exclude=cu01

module load miniconda/4.9.2
source activate smetana_3.7

id=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../matrix.txt)
mkdir -p /home/dongbiao/word_embedding_microbiome/modelseed/smetana/input/run_${id}
inputpath=/home/dongbiao/word_embedding_microbiome/modelseed/smetana/input/run_${id}
python run_smetana.py ../otu_pairs/${id}.csv ${inputpath}
rm -rf /home/dongbiao/word_embedding_microbiome/modelseed/smetana/input/run_${id}
