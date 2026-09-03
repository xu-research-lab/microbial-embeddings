#!/bin/bash

#SBATCH --job-name=glove_source
#SBATCH -N 1  
#SBATCH -p cu  
#SBATCH -n 28                 # Use 28 CPU cores
#SBATCH --mem=250G  
#SBATCH -o master_%a_m.log#SBATCH -e master_%a_m.err#SBATCH --array=1-8%3
#SBATCH --exclude=cu01

# Activate the Conda environment
source /home/cjj/miniconda3/etc/profile.d/conda.sh
conda activate microbiome_deep || {
  echo "Could not activate the Conda environment"
  exit 1
}
which python


# Read the parameter matrix
mapfile -t params < matrix_newdata_source.txt
IFS=' ' read -r p num <<< "${params[$SLURM_ARRAY_TASK_ID-1]}"

# id=$(sed "${SLURM_ARRAY_TASK_ID}q;d" matrix.txt)

# Global configuration
GENERATE_COOC=1                  # 0=reuse existing files, 1=generate new files
SOURCE_DIR="p80_abundance_percentile_1"  # Source directory when reusing files
BASE_DIR="."           # Base output directory
INPUT_BIOM="../../../data/gut_pretraining.biom"
CPUS_PER_TASK=28                       # CPU cores per task
# Create the base directory
mkdir -p "${BASE_DIR}"



# Process one percentile
process_percentile() {
  # TODO Fixed distance parameter: 80
  local percentile=$1
  local run=$2
  local work_dir="${BASE_DIR}/p${percentile}_${run}"
  local result_dir="${work_dir}/result"
  local log_dir="${work_dir}/log"
  # TODO: distance metric parameter
  local metric=${run}

  rm -rf "${result_dir}"
  rm -rf "${log_dir}"
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

    # Step 1: Generate the feature dictionary
    echo "Step 1/3: Generate the feature dictionary" | tee -a "${log_dir}/process.log"
    membed dict -b "${work_dir}/input.biom" -d "${work_dir}/feature-dict.csv" \
      > "${log_dir}/step1.log" 2>&1 || {
      echo "Feature dictionary generation failed. Percentile: ${percentile}" | tee -a "${log_dir}/error.log"
      return 1
    }

    # Step 2: Generate the co-occurrence matrix with --percentile_num
    echo "Step 2/3: Generate the co-occurrence matrix (percentile_num=${percentile})" | tee -a "${log_dir}/process.log"
    HDF5_USE_FILE_LOCKING=FALSE membed cooccur \
      -b "${work_dir}/input.biom" \
      -c "${work_dir}/table.co" \
      --metric ${metric} \
      --cpus ${CPUS_PER_TASK} \
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
  membed build-x-max-file -c "${work_dir}/table.co" -x "${work_dir}/xmax_file.npy" \
    --percentile_num ${percentile} \
    > "${log_dir}/step3b.log" 2>&1 || {
    echo "x-max file generation failed. Percentile: ${percentile}" | tee -a "${log_dir}/error.log"
    return 1
  }
  # Step 4: Train the GloVe model
  echo "Step 3/3: Train the model" | tee -a "${log_dir}/process.log"
  export OMP_NUM_THUMBREADS=${CPUS_PER_TASK}
  membed glove-train -d "${work_dir}/feature-dict.csv" \
    -c "${work_dir}/table.co" \
    -r "${result_dir}" \
    -x "${work_dir}/xmax_file.npy" \
    --lr 0.05 \
    --embedding-size 100 \
    --iter 100 \
    --cpus ${CPUS_PER_TASK} \
    > "${log_dir}/step3.log" 2>&1 || {
    echo "Model training failed. Percentile: ${percentile}" | tee -a "${log_dir}/error.log"
    return 1
  }

  echo "[$(date)] Percentile ${percentile}% completed" | tee -a "${log_dir}/process.log"
  return 0
}

# Process all percentiles sequentially
# exit_status=0
# for p in "${PERCENTILES[@]}"; do
#   for run in "${REPEATEDNUM[@]}"; do  # Run each percentile n times
#     if ! process_percentile ${p} ${run}; then
#       echo "[Error] ${p}% run ${run} failed"
#       exit_status=1
#     fi
#   done
# done

# Run one array task
if process_percentile $p $num; then
  echo "Task succeeded: p${p}_${num}"
else
  echo "Task failed: p${p}_${num}" >&2
  exit 1
fi

# Final status check
# if [ "${exit_status}" -eq 0 ]; then
#   echo "All percentile runs completed successfully"
# else
#   echo "Some runs failed; check the logs:"
#   find "${BASE_DIR}" -name "error.log" -exec grep -l "ERROR" {} \;
# fi

# exit ${exit_status}