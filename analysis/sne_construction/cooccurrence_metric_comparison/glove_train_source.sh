#!/bin/bash

#SBATCH --job-name=glove_array
#SBATCH -N 1
#SBATCH -p cu
#SBATCH -n 28                  # Use 28 CPU cores
#SBATCH --mem=100G
#SBATCH -o log/%x_%a.out       # Slurm standard output log
#SBATCH -e log/%x_%a.err       # Slurm error log
#SBATCH --array=1-35%3           # Submit 40 tasks (5 datasets * 8 metrics)
#SBATCH --exclude=cu01

# --- Global configuration ---
BASE_OUTPUT_DIR="./" # Root output directory for all results
PARAM_FILE="param_matrix.txt"    # Parameter matrix file
CPUS_PER_TASK=28                 # CPU cores per task
PERCENTILE_FOR_XMAX=80           # Percentile used by 'build-x-max-file'

# --- Activate environment ---
echo "Loading Conda environment..."
source /home/cjj/miniconda3/etc/profile.d/conda.sh
conda activate microbiome_deep || {
  echo "Error: Could not activate Conda environment 'microbiome_deep'"
  exit 1
}
echo "Python environment: $(which python)"
echo "==============================================="


# --- Read parameters for the current task ---
# Read row $SLURM_ARRAY_TASK_ID from the parameter matrix
mapfile -t params < ${PARAM_FILE}
IFS=' ' read -r biom_path metric <<< "${params[$SLURM_ARRAY_TASK_ID-1]}"

# Extract the dataset name from the BIOM path
dataset_name=$(basename "$biom_path" .biom)


# --- Set working and log directories ---
# Create a unique output directory for each metric and dataset
work_dir="${BASE_OUTPUT_DIR}/${metric}/${dataset_name}"
result_dir="${work_dir}/result_vectors" # Final embedding vectors
job_log_dir="${work_dir}/logs"          # Detailed task logs
local_biom_path="${work_dir}/input.biom"  # BIOM path in the working directory

# Recreate directories
rm -rf "${work_dir}"
mkdir -p "${result_dir}" "${job_log_dir}"


# --- Start task ---
echo "[$(date)] ==> Starting task #${SLURM_ARRAY_TASK_ID}"
echo "  - Input dataset: ${biom_path}"
echo "  - Metric      : ${metric}"
echo "  - Working directory: ${work_dir}"
echo "==============================================="

# Copy the BIOM file to the working directory
echo "Step 0/5: Copy the BIOM file to the working directory..." | tee -a "${job_log_dir}/process.log"
if [[ ! -f "${biom_path}" ]]; then
  echo "Error: Input file does not exist: ${biom_path}" | tee -a "${job_log_dir}/error.log"
  exit 1
fi
cp "${biom_path}" "${local_biom_path}" || {
    echo "Error: Failed to copy the BIOM file." | tee -a "${job_log_dir}/error.log"
    exit 1
}

# Step 1: Generate the feature dictionary
echo "Step 1/5: Generate the feature dictionary..." | tee -a "${job_log_dir}/process.log"
membed dict -b "${local_biom_path}" -d "${work_dir}/feature-dict.csv" \
  > "${job_log_dir}/step1_dict.log" 2>&1 || {
  echo "Error: Feature dictionary generation failed." | tee -a "${job_log_dir}/error.log"
  exit 1
}

# Step 2: Generate the co-occurrence matrix
echo "Step 2/5: Generate the co-occurrence matrix (Metric: ${metric})..." | tee -a "${job_log_dir}/process.log"
HDF5_USE_FILE_LOCKING=FALSE membed cooccur \
  -b "${local_biom_path}" \
  -c "${work_dir}/table.co" \
  --metric "${metric}" \
  --cpus "${CPUS_PER_TASK}" \
  > "${job_log_dir}/step2_cooccur.log" 2>&1 || {
  echo "Error: Co-occurrence matrix generation failed." | tee -a "${job_log_dir}/error.log"
  exit 1
}

# Step 3: Generate the x-max file
echo "Step 3/5: Generate the x-max file (Percentile: ${PERCENTILE_FOR_XMAX}%)..." | tee -a "${job_log_dir}/process.log"
membed build-x-max-file -c "${work_dir}/table.co" -x "${work_dir}/xmax_file.npy" \
  --percentile_num ${PERCENTILE_FOR_XMAX} \
  > "${job_log_dir}/step3_xmax.log" 2>&1 || {
  echo "Error: x-max file generation failed." | tee -a "${job_log_dir}/error.log"
  exit 1
}

# Step 4: Train the GloVe model
echo "Step 4/5: Train the GloVe model..." | tee -a "${job_log_dir}/process.log"
export OMP_NUM_THREADS=${CPUS_PER_TASK}
membed glove-train -d "${work_dir}/feature-dict.csv" \
  -c "${work_dir}/table.co" \
  -r "${result_dir}" \
  -x "${work_dir}/xmax_file.npy" \
  --lr 0.05 \
  --embedding-size 100 \
  --iter 100 \
  --cpus "${CPUS_PER_TASK}" \
  > "${job_log_dir}/step4_glove_train.log" 2>&1 || {
  echo "Error: Model training failed." | tee -a "${job_log_dir}/error.log"
  exit 1
}

# Step 5: Finish
echo "Step 5/5: Clean up and finish" | tee -a "${job_log_dir}/process.log"
echo "==============================================="
echo "[$(date)] ==> Task #${SLURM_ARRAY_TASK_ID} completed"
echo "Results saved to: ${result_dir}"
echo "==============================================="

exit 0