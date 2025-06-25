#!/bin/bash
#SBATCH --job-name=carvem
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -p cu
#SBATCH --mem=250
#SBATCH -o ../output/carvem_%a.out
#SBATCH -e ../output/carvem_%a.err
#SBATCH --array=1-10
#SBATCH --exclude=cu01,cu02,cu03,cu04,cu05,cu06,cu07,cu08

module load miniconda/4.9.2
source activate carvem_3.7

id=$(sed "${SLURM_ARRAY_TASK_ID}q;d" split_files.txt)
for i in `less split_files_add/${id}`;do
    python build_metbolic_model_picrust2.py ${i}
done
