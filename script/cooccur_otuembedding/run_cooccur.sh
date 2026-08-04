#!/usr/bin/env bash
#SBATCH --job-name=cooccur
#SBATCH -N 1  
#SBATCH -p cu  
#SBATCH -n 28                 
#SBATCH --mem=250G  
#SBATCH -o cooccur_%a.log  
#SBATCH -e cooccur_%a.err  
#SBATCH --array=1-8%4
#SBATCH --exclude=cu01

conda activate membed
data_path=../../data

# get the dict of OTU freature ID
# membed dict -b ${data_path}/gut_pretraining.biom -d feature-dict.csv

metric=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metrix.txt | tr -d '\r\n')
# caculate cooccur
mkdir ${metric}
membed cooccur -b ${data_path}/gut_pretraining.biom -c ${metric}/table.co --metric ${metric} --cpus 28