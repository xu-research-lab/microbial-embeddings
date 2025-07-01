#!/bin/bash  

#SBATCH --job-name=attention  
#SBATCH -N 1  
#SBATCH -p gpu  
#SBATCH -w gpu01  
#SBATCH -n 20       # 使用20核
#SBATCH --mem=100G  
#SBATCH -o log/4_21_attention_IBD_CRC.out  
#SBATCH -e log/4_21_attention_IBD_CRC.err  

source /home/cjj/miniconda3/etc/profile.d/conda.sh
conda activate jupyter_microbiome_deep

# 全局变量  
num_steps=600  
p_drop=0.4  
d_ff=8  
batch_size=512
# batch_size=32
d_model=100  
n_layers=1  
n_heads=1  
lr=0.001  
weight_decay=0.0001  
epoch=100
glove_embedding="transfer_p80_abundance_percentile_best_100_iter"
base_dir="/home/cjj/projects/memebed_analysis/result_new"
base_result_dir="./result"

# study_file="/home/dongbiao/all_study/data/existing_studies_all.txt"  
# mapfile -t studies < "$study_file" 

# 创建并行控制管道  
fifo="/tmp/$$.fifo"  
mkfifo "$fifo"  
exec 9<>"$fifo"  
rm -f "$fifo"  

# 初始化40个任务槽位（8GPU × 5任务/GPU）  
for gpu in {0..7}; do  
    for _ in {1..1}; do  
        echo "$gpu" >&9  
    done  
done  

for number in {16..20}; do  
    for disease_dir in "$base_dir"/*; do

        if [ ! -d "$disease_dir" ]; then
            continue
        fi
        disease=$(basename "$disease_dir")
        echo "Processing disease: $disease"

        # 跳过非IBD或CRC疾病
        # if [ "$disease" != "IBD_feces" ] && [ "$disease" != "CRC" ]; then
        #     # 这行代码块只会在 disease 既不是 IBD 也不是 CRC 的时候执行
        #     echo "Skipping $disease as it is not IBD or CRC." # 示例输出
        #     continue
        # fi
        if [ "$disease" != "IBD_feces" ] && [ "$disease" != "CRC" ]; then
            # 这行代码块只会在 disease 既不是 IBD 也不是 CRC 的时候执行
            echo "Skipping $disease as it is not IBD or CRC." # 示例输出
            continue
        fi
        disease="${disease}_${number}"

        # if [ "$disease" != "IBD" ] ; then
        #     # 这行代码块只会在 disease 既不是 IBD 也不是 CRC 的时候执行
        #     echo "Skipping $disease as it is not IBD or CRC." # 示例输出
        #     continue
        # fi



        # 获取该疾病下所有研究ID
        studies=()
        for study_dir in "$disease_dir"/*; do
            if [ -d "$study_dir" ]; then
                studies+=("$(basename "$study_dir")")
            fi
        done

        # 跳过只有一个研究的疾病
        if [ ${#studies[@]} -le 1 ]; then
            echo "Skipping $disease (only ${#studies[@]} studies)"
            continue
        fi

        # 遍历每个研究进行LOO训练
        for study_id in "${studies[@]}"; do
            # 获取GPU编号
            read -u9 numb

            (
                echo "[${disease}/${study_id}] => GPU-${numb}"
                data_dir="$disease_dir/$study_id"
                result_dir="$base_result_dir/$disease"
                mkdir -p "$result_dir"
                result_dir="$result_dir/$study_id/results"  # 结果存储在study_id子目录下的results文件夹
                
                # 创建结果目录
                mkdir -p "$result_dir"

                membed class-attention -g "suffle_embedding/${glove_embedding}_${number}" \
                    -tra_otu "$data_dir/train_loo.biom" \
                    -tes_otu "$data_dir/test_loo.biom" \
                    -m "$disease_dir/metadata.tsv" \
                    -ploss "${result_dir}/loss_loo.png" \
                    -pauc "${result_dir}/auc_loo.png" \
                    -e "${result_dir}/attention_loo.pt" \
                    --num-steps $num_steps \
                    --group group \
                    --loss BCE_loss \
                    --p-drop $p_drop \
                    --d-ff $d_ff \
                    --batch-size $batch_size \
                    --d-model $d_model \
                    --n-layers $n_layers \
                    --n-heads $n_heads \
                    --numb $numb \
                    --lr $lr \
                    --weight-decay $weight_decay \
                    --num-epochs $epoch

                # 释放GPU槽位
                echo "$numb" >&9
            ) &
        done
    done 
done



# 等待所有后台任务完成  
wait  
exec 9>&-  