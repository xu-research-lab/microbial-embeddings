#!/usr/bin/env bash

#SBATCH --job-name=glove
#SBATCH -N 1  
#SBATCH -p cu  
#SBATCH -n 28                 
#SBATCH --mem=250G  
#SBATCH -o glove_%a.log  
#SBATCH -e glove_%a.err  
#SBATCH --array=3
#SBATCH --exclude=cu01

module load miniconda/4.9.2 
source activate membed

metric=$(sed -n "${SLURM_ARRAY_TASK_ID}p" metrix.txt | tr -d '\r\n')
# Pick the x_max weighting cutoff from the co-occurrence records
membed build-x-max-file -c ${metric}/table.co -x ${metric}/xmax_file.npy --percentile-num 80

# Run GloVe embedding
membed glove-train \
       -d feature-dict.csv \
       -c ${metric}/table.co \
       -r ${metric} \
       -x ${metric}/xmax_file.npy \
       --lr 0.05 \
       --embedding-size 100 \
       --iter 100 --cpus 28