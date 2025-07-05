#!/bin/bash
#SBATCH --job-name=smetana
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p cu
#SBATCH --mem=100
#SBATCH -o smetana_%a.out
#SBATCH -e smetana_%a.err
#SBATCH --array=1
#SBATCH --exclude=cu01

module load miniconda/4.9.2
source activate smetana_3.7

id=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../matrix.txt)
mkdir -p Data/smetana/input/run_${id}
inputpath=Data/smetana/input/run_${id}
python run_smetana.py otu_pairs/${id}.csv ${inputpath}
rm -rf Data/smetana/input/run_${id}
