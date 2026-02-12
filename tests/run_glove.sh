#!/bin/bash

# Script to run Glove pipeline on test_Glove.biom
# This script follows the steps outlined in the README

set -e  # Exit on error

# Configuration
INPUT_BIOM="data/test_Glove.biom"
OUTPUT_DIR="glove_output"
METRIC="abundance_percentile"  # Metric used in the study
PERCENTILE_FOR_XMAX=80
CPUS=8  # Adjust based on your system
EMBEDDING_SIZE=100
ITERATIONS=100
LEARNING_RATE=0.05
LOG_FILE="${OUTPUT_DIR}/pipeline_timing.log"

# Time logging function
log_time() {
    local step="$1"
    local status="$2"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] ${step}: ${status}" | tee -a "${LOG_FILE}"
}

# Elapsed time calculation function
print_elapsed() {
    local step="$1"
    local start_time="$2"
    local end_time="$3"
    local elapsed=$((end_time - start_time))
    local hours=$((elapsed / 3600))
    local minutes=$(( (elapsed % 3600) / 60 ))
    local seconds=$((elapsed % 60))
    printf "[%s] %s completed in: %02d:%02d:%02d (hh:mm:ss)\n" \
        "$(date '+%Y-%m-%d %H:%M:%S')" \
        "$step" "$hours" "$minutes" "$seconds" | tee -a "${LOG_FILE}"
}

echo "==============================================="
echo "Running Glove pipeline on test_Glove.biom"
echo "==============================================="

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Initialize log file
echo "===============================================" > "${LOG_FILE}"
echo "Glove Pipeline Timing Log - $(date '+%Y-%m-%d %H:%M:%S')" >> "${LOG_FILE}"
echo "===============================================" >> "${LOG_FILE}"
echo "Configuration:" >> "${LOG_FILE}"
echo "  - CPUs: ${CPUS}" >> "${LOG_FILE}"
echo "  - Embedding size: ${EMBEDDING_SIZE}" >> "${LOG_FILE}"
echo "  - Iterations: ${ITERATIONS}" >> "${LOG_FILE}"
echo "===============================================" >> "${LOG_FILE}"

# Check if input file exists
if [[ ! -f "${INPUT_BIOM}" ]]; then
    echo "Error: Input file not found: ${INPUT_BIOM}"
    exit 1
fi

# Step 1: Generate feature dictionary
log_time "Step 1/5" "Started - Generating feature dictionary"
step1_start=$(date +%s)
membed dict -b "${INPUT_BIOM}" -d "${OUTPUT_DIR}/feature-dict.csv"
step1_end=$(date +%s)
print_elapsed "Step 1 - Feature dictionary" $step1_start $step1_end
log_time "Step 1/5" "Completed - feature-dict.csv saved"

echo ""

# Step 2: Generate co-occurrence matrix
log_time "Step 2/5" "Started - Generating co-occurrence matrix (metric: ${METRIC})"
step2_start=$(date +%s)
HDF5_USE_FILE_LOCKING=FALSE membed cooccur \
    -b "${INPUT_BIOM}" \
    -c "${OUTPUT_DIR}/table.co" \
    --metric "${METRIC}" \
    --cpus "${CPUS}"
step2_end=$(date +%s)
print_elapsed "Step 2 - Co-occurrence matrix" $step2_start $step2_end
log_time "Step 2/5" "Completed - table.co saved"

echo ""

# Step 3: Generate xmax file
log_time "Step 3/5" "Started - Generating xmax file (percentile: ${PERCENTILE_FOR_XMAX}%)"
step3_start=$(date +%s)
membed build-x-max-file \
    -c "${OUTPUT_DIR}/table.co" \
    -x "${OUTPUT_DIR}/xmax_file.npy" \
    --percentile_num "${PERCENTILE_FOR_XMAX}"
step3_end=$(date +%s)
print_elapsed "Step 3 - xmax file" $step3_start $step3_end
log_time "Step 3/5" "Completed - xmax_file.npy saved"

echo ""

# Step 4: Train Glove model
log_time "Step 4/5" "Started - Training Glove model"
step4_start=$(date +%s)
export OMP_NUM_THREADS=${CPUS}
membed glove-train \
    -d "${OUTPUT_DIR}/feature-dict.csv" \
    -c "${OUTPUT_DIR}/table.co" \
    -r "${OUTPUT_DIR}" \
    -x "${OUTPUT_DIR}/xmax_file.npy" \
    --lr "${LEARNING_RATE}" \
    --embedding-size "${EMBEDDING_SIZE}" \
    --iter "${ITERATIONS}" \
    --cpus "${CPUS}"
step4_end=$(date +%s)
print_elapsed "Step 4 - Glove training" $step4_start $step4_end
log_time "Step 4/5" "Completed - embeddings_${EMBEDDING_SIZE}.txt saved"

echo ""

# Step 5: Complete
step5_start=$(date +%s)
log_time "Step 5/5" "Pipeline completed"
step5_end=$(date +%s)
print_elapsed "Total pipeline" $step1_start $step5_end

echo "==============================================="
echo "Embeddings saved to: ${OUTPUT_DIR}/embeddings_${EMBEDDING_SIZE}.txt"
echo "Timing log saved to: ${LOG_FILE}"
echo "==============================================="

# Display timing summary
echo ""
echo "==============================================="
echo "TIMING SUMMARY"
echo "==============================================="
grep "completed in:" "${LOG_FILE}"
echo "==============================================="