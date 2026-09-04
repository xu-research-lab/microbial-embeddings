#!/bin/bash
#SBATCH --job-name=glove_embedding_size
#SBATCH -N 1  
#SBATCH -p cu  
#SBATCH -n 28                 # Use 28 CPU cores
#SBATCH --mem=250G  
#SBATCH --exclude=cu01


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

# Activate the Conda environment
source /home/cjj/miniconda3/etc/profile.d/conda.sh
conda activate microbiome_deep || {
  echo "Could not activate the Conda environment"
  exit 1
}
which python

set -euo pipefail


# Read percentile, metric, and embedding size from command-line arguments.
p="${1:-80}"
num="${2:-abundance_percentile}"
embedding_size="${3:-100}"


# Global configuration
GENERATE_COOC="${GENERATE_COOC:-1}"  # 0=reuse existing files, 1=generate new files
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}/results}"
SOURCE_DIR="${SOURCE_DIR:-${BASE_DIR}/p80_${num}}"
INPUT_BIOM="${INPUT_BIOM:-${REPO_ROOT}/data/gut_pretraining.biom}"
CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-28}"



# Process one percentile
process_percentile() {
  local percentile=$1
  local run=$2
  local embedding_size="$3"
  local work_dir="${BASE_DIR}/p${percentile}_${run}_emb${embedding_size}"
  local result_dir="${work_dir}/result"
  local log_dir="${work_dir}/log"

  local metric="${run}"

  rm -rf "${result_dir}" "${log_dir}"
  mkdir -p "${work_dir}" "${result_dir}" "${log_dir}"

  echo "==============================================="
  echo "[$(date)] Starting percentile ${percentile}%"
  echo "Working directory: ${work_dir}"
  echo "Log directory: ${log_dir}"
  
  # Create the working directory and copy the BIOM file
  if [[ ! -f "${INPUT_BIOM}" ]]; then
    echo "Error: Input file does not exist: ${INPUT_BIOM}" | tee -a "${log_dir}/error.log"
    return 1
  fi
  cp "${INPUT_BIOM}" "${work_dir}/input.biom" || return 1

  # 0=reuse existing files, 1=generate a new co-occurrence matrix
  if [[ "${GENERATE_COOC}" -eq 1 ]]; then
    rm -f "${work_dir}/feature-dict.csv" "${work_dir}/table.co" "${work_dir}/xmax_file.npy"

    # Step 1: Generate the feature dictionary
    echo "Step 1/4: Generate the feature dictionary" | tee -a "${log_dir}/process.log"
    membed dict -b "${work_dir}/input.biom" -d "${work_dir}/feature-dict.csv" \
      > "${log_dir}/step1.log" 2>&1 || {
      echo "Feature dictionary generation failed. Percentile: ${percentile}" | tee -a "${log_dir}/error.log"
      return 1
    }

    # Step 2: Generate the co-occurrence matrix
    echo "Step 2/4: Generate the co-occurrence matrix (metric=${metric})" | tee -a "${log_dir}/process.log"
    HDF5_USE_FILE_LOCKING=FALSE membed cooccur \
      -b "${work_dir}/input.biom" \
      -c "${work_dir}/table.co" \
      --metric "${metric}" \
      --cpus "${CPUS_PER_TASK}" \
      > "${log_dir}/step2.log" 2>&1 || {
      echo "Co-occurrence matrix generation failed. Percentile: ${percentile}" | tee -a "${log_dir}/error.log"
      return 1
    }
  else
    # Reuse mode
    echo "Reusing existing files from: ${SOURCE_DIR}" | tee -a "${log_dir}/process.log"
    local required_files=(
      "${SOURCE_DIR}/feature-dict.csv"
      "${SOURCE_DIR}/table.co" 
    )
    
    # Check that source files exist
    for f in "${required_files[@]}"; do
      if [[ ! -f "${f}" ]]; then
        echo "Missing required file: ${f}" | tee -a "${log_dir}/error.log"
        return 1
      fi
      cp "${f}" "${work_dir}/" || return 1
    done
  fi

  # Step 3: Generate the x-max value
  echo "Step 3/4: Generate the x-max file" | tee -a "${log_dir}/process.log"
  membed build-x-max-file -c "${work_dir}/table.co" -x "${work_dir}/xmax_file.npy" \
    --percentile_num "${percentile}" \
    > "${log_dir}/step3b.log" 2>&1 || {
    echo "x-max file generation failed. Percentile: ${percentile}" | tee -a "${log_dir}/error.log"
    return 1
  }
  # Step 4: Train the GloVe model
  echo "Step 4/4: Train the model" | tee -a "${log_dir}/process.log"
  export OMP_NUM_THREADS="${CPUS_PER_TASK}"
  membed glove-train -d "${work_dir}/feature-dict.csv" \
    -c "${work_dir}/table.co" \
    -r "${result_dir}" \
    -x "${work_dir}/xmax_file.npy" \
    --lr 0.05 \
    --embedding-size "${embedding_size}" \
    --iter 100 \
    --cpus "${CPUS_PER_TASK}" \
    > "${log_dir}/step3.log" 2>&1 || {
    echo "Model training failed. Percentile: ${percentile}" | tee -a "${log_dir}/error.log"
    return 1
  }

  echo "[$(date)] Percentile ${percentile}% completed" | tee -a "${log_dir}/process.log"
  return 0
}

# Run one array task
if process_percentile "${p}" "${num}" "${embedding_size}"; then
  echo "Task succeeded: p${p}_${num}_${embedding_size}"
else
  echo "Task failed: p${p}_${num}_${embedding_size}" >&2
  exit 1
fi
