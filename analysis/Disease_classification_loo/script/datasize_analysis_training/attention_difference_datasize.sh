#!/bin/bash

#SBATCH --job-name=datasize_attention
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH -w gpu01
#SBATCH -n 20
#SBATCH --mem=50G
#SBATCH -o log/datasize_attention_%j.out  # Use job ID for unique logs
#SBATCH -e log/datasize_attention_%j.err  # Use job ID for unique logs

# Create log directory if it doesn't exist
mkdir -p log

# module load miniconda/4.9.2
# source activate jupyter_notebook
conda init bash
conda activate jupyter_microbiome_deep

# --- Global Variables ---
num_steps=600
p_drop=0.4
d_ff=8
batch_size=512
d_model=100  # Assuming this is constant based on _100.txt examples and original script
n_layers=1
n_heads=1
lr=0.001
weight_decay=0.0001
epoch=100

base_data_dir="/home/cjj/projects/memebed_analysis/result_new" # Input data directory (train/test/metadata)
embeding_base_dir="/home/cjj/projects/memebed_analysis/script/script_analysis_datasize/embedding" # Base dir for datasize_result_* folders
result_base_dir="/home/cjj/projects/memebed_analysis/script/script_analysis_datasize/result"   # Base dir for output results

target_diseases=("IBD_feces" "CRC") # Add other diseases here if needed
# target_diseases=("CRC") # Add other diseases here if needed
# --- Parallel GPU Management ---
fifo="/tmp/gpu_fifo_$$.fifo"
mkfifo "$fifo"
exec 9<>"$fifo" # Open FIFO for reading and writing on FD 9
rm -f "$fifo" # Remove fifo file node, the pipe persists

# Initialize GPU task slots (Assuming 8 GPUs on the node gpu01)
for gpu_id in {0..7}; do
    echo "$gpu_id" >&9
done

echo "Starting script..."
echo "Embedding base directory: $embeding_base_dir"
echo "Result base directory: $result_base_dir"

# --- Main Loop ---
for disease in "${target_diseases[@]}"; do
    disease_data_dir="$base_data_dir/$disease"
    echo "Processing disease: $disease"

    if [ ! -d "$disease_data_dir" ]; then
        echo "Warning: Disease data directory $disease_data_dir not found! Skipping."
        continue
    fi
    # if [ $disease != "CRC" ]; then
    #     echo "Warning: CRC data is too large to process. Skipping."
    #     continue
    # fi

    studies=()
    for study_dir in "$disease_data_dir"/*; do
        if [ -d "$study_dir" ]; then
            studies+=("$(basename "$study_dir")")
        fi
    done

    if [ ${#studies[@]} -le 1 ]; then
        echo "Skipping $disease (only ${#studies[@]} studies found, need > 1 for LOO)"
        continue
    fi
    echo "Found studies for $disease: ${studies[*]}"

    # Loop through data size result directories (datasize_result_1, datasize_result_2, etc.)
    for data_size_dir in "$embeding_base_dir"/datasize_result_*; do
        if [ ! -d "$data_size_dir" ]; then
            # echo "Skipping non-directory item: $data_size_dir" # Optional debug
            continue
        fi
        # --- Get the name of the data size directory (e.g., "datasize_result_1") ---
        data_size_dirname=$(basename "$data_size_dir")
        echo "Processing data size directory: $data_size_dirname"

        # Loop through embedding files within the current data size directory
        for embedding_file in "$data_size_dir"/subset_table_*.txt; do
            if [ ! -f "$embedding_file" ]; then
                 continue # Skip if not a valid file
            fi
            # if [ "$embedding_file" != "$data_size_dir/subset_table_1k_100.txt" ] && [ "$embedding_file" != "$data_size_dir/subset_table_2k_100.txt" ] && [ "$embedding_file" != "$data_size_dir/subset_table_5k_100.txt" ] ; then
            #     continue # Skip if not a valid file
            # fi
            # if [ "$embedding_file" != "$data_size_dir/subset_table_1k_100.txt" ] && [ "$embedding_file" != "$data_size_dir/subset_table_2k_100.txt" ] ; then
            #     continue # Skip if not a valid file
            # fi

            temp_path_no_ext="${embedding_file%.txt}" # Removes .txt -> /path/to/base_name_100
            membed_g_arg="${temp_path_no_ext%_*}"     # Removes _100 -> /path/to/base_nam

            # --- Extract base embedding name (e.g., "subset_table_1w_100") ---
            embedding_name=$(basename "$embedding_file" .txt)
            # --- Extract the short embedding name (e.g., "subset_table_1w") ---
            # This removes the last underscore and following characters (assuming format like _100)
            embedding_name_short="${embedding_name%_*}"
            echo "  Processing embedding: $embedding_name (short name: $embedding_name_short)"


            # Iterate through each study for Leave-One-Out (LOO) training/testing
            for study_id in "${studies[@]}"; do
                read -u9 gpu_num # Acquire a GPU ID

                # Run the training in the background
                (
                    # --- CONSTRUCT THE NEW RESULT PATH ---
                    # Format: /result_base_dir / Disease / DataSizeDirName / EmbeddingShortName / StudyID / results
                    loo_result_dir="$result_base_dir/$disease/$data_size_dirname/$embedding_name_short/$study_id"

                    echo "    [${disease}/${data_size_dirname}/${embedding_name_short}/${study_id}] Starting job on GPU-${gpu_num}. Output: $loo_result_dir"

                    # Create the specific result directory for this run *before* running the command
                    mkdir -p "$loo_result_dir"

                    loo_data_dir="$disease_data_dir/$study_id" # Dir containing train_loo.biom, test_loo.biom

                    # Check if necessary input files exist
                    if [ ! -f "$loo_data_dir/train_loo.biom" ] || \
                       [ ! -f "$loo_data_dir/test_loo.biom" ] || \
                       [ ! -f "$disease_data_dir/metadata.tsv" ]; then
                        echo "Error: Missing input files for ${disease}/${data_size_dirname}/${embedding_name_short}/${study_id}. Skipping."
                        echo "$gpu_num" >&9 # Release GPU slot
                        exit 1 # Exit subshell
                    fi

                    # --- Run the membed command ---
                    membed class-attention -g "$membed_g_arg" \
                        -tra_otu "$loo_data_dir/train_loo.biom" \
                        -tes_otu "$loo_data_dir/test_loo.biom" \
                        -m "$disease_data_dir/metadata.tsv" \
                        -ploss "${loo_result_dir}/loss_loo.png" \
                        -pauc "${loo_result_dir}/auc_loo.png" \
                        -e "${loo_result_dir}/attention_loo.pt" \
                        --num-steps $num_steps \
                        --group group \
                        --loss BCE_loss \
                        --p-drop $p_drop \
                        --d-ff $d_ff \
                        --batch-size $batch_size \
                        --d-model $d_model \
                        --n-layers $n_layers \
                        --n-heads $n_heads \
                        --numb $gpu_num \
                        --lr $lr \
                        --weight-decay $weight_decay \
                        --num-epochs $epoch

                    exit_status=$?
                    if [ $exit_status -ne 0 ]; then
                       echo "    [${disease}/${data_size_dirname}/${embedding_name_short}/${study_id}] Job on GPU-${gpu_num} failed with status $exit_status."
                    else
                       echo "    [${disease}/${data_size_dirname}/${embedding_name_short}/${study_id}] Job on GPU-${gpu_num} completed."
                    fi

                    # Release the GPU slot
                    echo "$gpu_num" >&9

                ) & # Run the subshell in the background

            done # End study loop (LOO)
        done # End embedding file loop
    done # End data size dir loop
done # End disease loop

# Wait for all background jobs
wait
echo "All background jobs finished."

# Close the file descriptor
exec 9>&-
echo "Script finished."