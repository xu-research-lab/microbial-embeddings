#!/bin/bash

# Script to run attention-based classification on IBD data

set -e  # Exit on error

# Configuration
EMBEDDING_FILE="../data/social_niche_embedding_100.txt"
TRAIN_BIOM="data/IBD_train.biom"
VALID_BIOM="data/IBD_valid.biom"
TEST_BIOM="data/IBD_test.biom"
METADATA_FILE="data/metadata_IBD.txt"
OUTPUT_DIR="classification_output"
LOG_FILE="${OUTPUT_DIR}/classification_time.txt"

# Hyperparameters
LABELS_COL="group"
SAMPLE_ID_COL="sample"
NUM_STEPS=400
NUM_EPOCHS=100
LOSS="BCE_loss"
DROPOUT=0.4
D_FF=8
D_MODEL=100
N_LAYERS=1
N_HEADS=1
WEIGHT_DECAY=0.0001
LEARNING_RATE=0.001
BATCH_SIZE=64
# Zero-based CUDA device index passed to --numb; falls back to CPU when no
# CUDA device is available. Override with e.g. NUMB=-1 to force CPU.
NUMB="${NUMB:-0}"

echo "==============================================="
echo "Running attention-based classification on IBD data"
echo "==============================================="

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Check if required files exist
for file in "${EMBEDDING_FILE}" "${TRAIN_BIOM}" "${VALID_BIOM}" "${TEST_BIOM}" "${METADATA_FILE}"; do
    if [[ ! -f "${file}" ]]; then
        echo "Error: Required file not found: ${file}"
        exit 1
    fi
done

echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')" | tee "${LOG_FILE}"
echo "Command: membed class-attention ..." >> "${LOG_FILE}"
echo "" >> "${LOG_FILE}"

# Start timer
start_time=$(date +%s)

# Run the classification command
HDF5_USE_FILE_LOCKING=FALSE membed class-attention \
    -g "${EMBEDDING_FILE}" \
    -tra_otu "${TRAIN_BIOM}" \
    -val_otu "${VALID_BIOM}" \
    -tes_otu "${TEST_BIOM}" \
    -m "${METADATA_FILE}" \
    --labels_col "${LABELS_COL}" \
    --sample_id_col "${SAMPLE_ID_COL}" \
    -e "${OUTPUT_DIR}/attention_weights.pt" \
    -ploss "${OUTPUT_DIR}/attention_loss.png" \
    -pauc "${OUTPUT_DIR}/attention_auc.png" \
    --num-steps "${NUM_STEPS}" \
    --num-epochs "${NUM_EPOCHS}" \
    --loss "${LOSS}" \
    --p-drop "${DROPOUT}" \
    --d-ff "${D_FF}" \
    --d-model "${D_MODEL}" \
    --n-layers "${N_LAYERS}" \
    --n-heads "${N_HEADS}" \
    --weight-decay "${WEIGHT_DECAY}" \
    --lr "${LEARNING_RATE}" \
    --batch-size "${BATCH_SIZE}" \
    --numb "${NUMB}"

# End timer
end_time=$(date +%s)
elapsed=$((end_time - start_time))
hours=$((elapsed / 3600))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$((elapsed % 60))

echo "" | tee -a "${LOG_FILE}"
echo "End time: $(date '+%Y-%m-%d %H:%M:%S')" | tee -a "${LOG_FILE}"
echo "----------------------------------------" | tee -a "${LOG_FILE}"
printf "Total execution time: %02d:%02d:%02d (hh:mm:ss)\n" $hours $minutes $seconds | tee -a "${LOG_FILE}"
printf "Total seconds: %d seconds\n" $elapsed | tee -a "${LOG_FILE}"
echo "----------------------------------------" | tee -a "${LOG_FILE}"
echo "Output saved to: ${OUTPUT_DIR}" | tee -a "${LOG_FILE}"
echo "Time log saved to: ${LOG_FILE}" | tee -a "${LOG_FILE}"
echo "==============================================="