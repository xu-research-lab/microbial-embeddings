#!/bin/bash
#SBATCH --job-name=diamond
#SBATCH -N 1
#SBATCH -n 5
#SBATCH -p cu
#SBATCH --mem=40
#SBATCH -o diamond_%a.out
#SBATCH -e diamond_%a.err
#SBATCH --array=2-100
#SBATCH --exclude=cu01

module load miniconda/4.9.2
source activate carvem_3.7

id=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../split_file_genome)
for i in `less ../split_files/${id}`;do
    python run_diamond.py ${i}
done