#!/bin/bash
#SBATCH --job-name=rf
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH -w gpu01
#SBATCH -n 4
#SBATCH --mem=50G
#SBATCH -o log/4_21_RF.out
#SBATCH -e log/4_21_RF.err

source /home/cjj/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_microbiome_deep

script_path="." 
base_dir="./data"

# 遍历疾病目录
for disease_dir in "$base_dir"/*; do
    # 跳过非目录文件
    [ -d "$disease_dir" ] || continue
    
    disease=$(basename "$disease_dir")
    echo "Processing disease: $disease"
    
    # 获取研究列表
    studies=()
    while IFS= read -d $'\0' -r study; do
        studies+=("$(basename "$study")")
    done < <(find "$disease_dir" -mindepth 1 -maxdepth 1 -type d -print0)

    # 检查研究数量
    if [ ${#studies[@]} -le 1 ]; then
        echo "Skipping $disease (only ${#studies[@]} studies)"
        continue
    fi

    # 处理每个研究
    for study_id in "${studies[@]}"; do
        data_dir="$disease_dir/$study_id"
        result_dir="$data_dir/RF"
        
        # 创建结果目录
        if ! mkdir -p "$result_dir"; then
            echo "Failed to create directory: $result_dir"
            continue
        fi

        # 检查数据文件是否存在
        if [ ! -f "$data_dir/train_loo.biom" ] || [ ! -f "$data_dir/test_loo.biom" ]; then
            echo "Missing data files in $data_dir"
            continue
        fi

        echo "Processing $disease/$study_id"

        echo ${study_id} K_ >&2
        echo ${study_id} K_
        # 运行分类器
        python "Randomforestclassifier.py" \
            "$data_dir/train_loo.biom" \
            "$data_dir/test_loo.biom" \
            "$disease_dir/metadata.tsv" \
            None \
            sex \
            "$result_dir/loo.png" \
            "$result_dir/RF_ROC.csv" \
            "$result_dir/RF_Scores.csv" \
            "$result_dir/Feature_importance.csv"
    done
done