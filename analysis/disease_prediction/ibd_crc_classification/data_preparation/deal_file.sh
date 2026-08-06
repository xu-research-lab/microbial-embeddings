#!/bin/bash

# --- CONFIGURATION ---
# The top-level directory containing the study folders (e.g., PRJNA324147)
SOURCE_DIR="/home/cjj/projects/memebed_analysis/result_new/CRC"

# The destination folder for all the renamed files
DEST_DIR="/home/cjj/projects/memebed_analysis/script/script_leave_one_out_new/CRC_feces_attention"

# --- SCRIPT LOGIC ---

# 1. Ensure the destination directory exists
mkdir -p "$DEST_DIR"

# 2. Initialize a counter for sequential renaming
i=1

echo "Starting file copy and rename process..."

# 3. Loop through each study directory within the source directory
# The trailing slash on '*/' ensures we only loop through directories
for study_path in "$SOURCE_DIR"/*/; do

    # Define the full paths to the three files we need to copy
    attention_src="${study_path}results/attention_loo.pt"
    test_src="${study_path}test_loo.biom"
    train_src="${study_path}train_loo.biom"

    # Check if all three source files exist before proceeding
    if [ -f "$attention_src" ] && [ -f "$test_src" ] && [ -f "$train_src" ]; then
    
        echo "Copying set $i from $(basename "$study_path")"

        # 4. Copy each file to the destination, renaming it with the counter 'i'
        cp "$attention_src" "${DEST_DIR}/attention_${i}.pt"
        cp "$test_src"      "${DEST_DIR}/test_${i}.biom"
        cp "$train_src"     "${DEST_DIR}/train_${i}.biom"

        # 5. Increment the counter for the next study
        i=$((i + 1))
        
    else
        # Print a warning if a study directory is missing one or more files
        echo "⚠️  Skipping $(basename "$study_path") - required files not found."
    fi
done

echo ""
echo "✅  Process complete. Total sets of files copied: $((i - 1))."
echo "Files are located in: $DEST_DIR"