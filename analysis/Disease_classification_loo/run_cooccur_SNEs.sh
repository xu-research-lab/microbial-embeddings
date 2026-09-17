#!/bin/bash
#SBATCH --job-name=cooccur
#SBATCH -N 1
#SBATCH -p cu
#SBATCH -n 28
#SBATCH --mem=250G
#SBATCH -o cooccur_%a.log
#SBATCH -e cooccur_%a.err
#SBATCH --array=1-5
#SBATCH --exclude=cu01,cu02

# Pre-training data-size sweep: one GloVe run per (replicate, subset size).
# Inputs come from get_subdatasets.py; the array index is the replicate.

set -euo pipefail

source /home/dongbiao/miniconda3/etc/profile.d/conda.sh
conda activate membed

rep=${SLURM_ARRAY_TASK_ID}

data_path=Data/pretraining_datasize/trainning_data/data_${rep}
output_path=Data/pretraining_datasize/embedding/datasize_result_${rep}
mkdir -p "${output_path}"

# Same sizes get_subdatasets.py writes.
SUBSETS=("1k" "2k" "5k" "1w" "2w" "4w" "8w" "16w")

EMB_SIZE=100
ITER=100
LR=0.05
CPUS=28
PERCENTILE_NUM=80

for s in "${SUBSETS[@]}"; do
    echo "=== Processing subset ${s} ==="

    biom_file="${data_path}/subset_table_${s}.biom"

    # get the dict of OTU feature ID
    membed dict -b "${biom_file}" -d "${output_path}/feature-dict_${s}.csv"

    # cooccurrence
    membed cooccur -b "${biom_file}" -c "${output_path}/table_${s}.co" \
           --metric abundance_percentile --cpus "${CPUS}"

    # pick the x_max weighting cutoff from the co-occurrence records
    membed build-x-max-file -c "${output_path}/table_${s}.co" \
           -x "${output_path}/xmax_file_${s}.npy" --percentile-num "${PERCENTILE_NUM}"

    # run GloVe embedding
    membed glove-train \
           -d "${output_path}/feature-dict_${s}.csv" \
           -c "${output_path}/table_${s}.co" \
           -r "${output_path}/subset_table_${s}" \
           -x "${output_path}/xmax_file_${s}.npy" \
           --lr "${LR}" \
           --embedding-size "${EMB_SIZE}" \
           --iter "${ITER}" --cpus "${CPUS}"

    # the file the downstream classifier reads (and the one tracked in git)
    cp "${output_path}/subset_table_${s}/embeddings_${EMB_SIZE}.txt" \
       "${output_path}/subset_table_${s}_${EMB_SIZE}.txt"
done
