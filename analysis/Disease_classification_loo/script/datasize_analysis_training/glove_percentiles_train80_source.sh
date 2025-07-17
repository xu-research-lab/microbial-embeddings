#!/bin/bash

#SBATCH --job-name=glove_source
#SBATCH -N 1  
#SBATCH -p cu  
#SBATCH -n 28                 # 使用28核
#SBATCH --mem=100G  
#SBATCH -o log/master_%a_r.log  # 改动点
#SBATCH -e log/master_%a_r.err  # 改动点
#SBATCH --array=1-3
#SBATCH --exclude=cu01,

# 激活conda环境
source /home/cjj/miniconda3/etc/profile.d/conda.sh
conda activate microbiome_deep || {
  echo "无法激活conda环境"
  exit 1
}
which python


# 读取参数矩阵
mapfile -t params < dataset_list.txt # 改动点
IFS=' ' read -r path <<< "${params[$SLURM_ARRAY_TASK_ID-1]}"

# id=$(sed "${SLURM_ARRAY_TASK_ID}q;d" matrix.txt)

# 全局配置
CPUS_PER_TASK=28                       # 每个任务CPU数

# 处理单个百分位数的函数
process_percentile() {
  # TODO 固定距离调参 80
  local percentile=80
  local bim_path=$1
  local parent_dir=$(dirname "$bim_path")
  local work_dir="${parent_dir}" # 改动点
  local result_dir="${work_dir}/result"
  local log_dir="${work_dir}/log"
  # TODO 距离调参
  local metric="abundance_percentile"

  rm -rf "${result_dir}"
  rm -rf "${log_dir}"
  mkdir -p "${work_dir}" "${result_dir}" "${log_dir}"

  echo "==============================================="
  echo "[$(date)] 开始处理百分位数 ${percentile}%"
  echo "工作目录: ${work_dir}"
  echo "日志目录: ${log_dir}"
  
  # 创建工作目录并复制biom文件
  # if [[ ! -f "${INPUT_BIOM}" ]]; then
  #   echo "错误：输入文件不存在 ${INPUT_BIOM}" | tee -a "${log_dir}/error.log"
  #   return 1
  # fi
  # cp "${INPUT_BIOM}" "${work_dir}/input.biom" || return 1

  # Step 1: 生成特征字典
  echo "步骤1/3：生成特征字典" | tee -a "${log_dir}/process.log"
  membed dict -b "${work_dir}/input.biom" -d "${work_dir}/feature-dict.csv" \
    > "${log_dir}/step1.log" 2>&1 || {
    echo "字典生成失败！百分位数: ${percentile}" | tee -a "${log_dir}/error.log"
    return 1
  }

  # Step 2: 生成共现矩阵（关键修改：添加--percentile_num参数）
  echo "步骤2/3：生成共现矩阵 (percentile_num=${percentile})" | tee -a "${log_dir}/process.log"
  HDF5_USE_FILE_LOCKING=FALSE membed cooccur \
    -b "${work_dir}/input.biom" \
    -c "${work_dir}/table.co" \
    --metric ${metric} \
    --cpus ${CPUS_PER_TASK} \
    > "${log_dir}/step2.log" 2>&1 || {
    echo "共现矩阵生成失败！百分位数: ${percentile}" | tee -a "${log_dir}/error.log"
    return 1
  }

  # Step 3: 生成xmax值
  membed build-x-max-file -c "${work_dir}/table.co" -x "${work_dir}/xmax_file.npy" \
    --percentile_num ${percentile} \
    > "${log_dir}/step3b.log" 2>&1 || {
    echo "xmax文件生成失败！百分位数: ${percentile}" | tee -a "${log_dir}/error.log"
    return 1
  }
  # Step 4: 训练Glove模型
  echo "步骤3/3：训练模型" | tee -a "${log_dir}/process.log"
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
    echo "模型训练失败！百分位数: ${percentile}" | tee -a "${log_dir}/error.log"
    return 1
  }

  echo "[$(date)] 百分位数 ${percentile}% 处理完成" | tee -a "${log_dir}/process.log"
  return 0
}

# 顺序执行所有百分位数处理
# exit_status=0
# for p in "${PERCENTILES[@]}"; do
#   for run in "${REPEATEDNUM[@]}"; do  # 每个百分位数运行n次
#     if ! process_percentile ${p} ${run}; then
#       echo "[错误] ${p}% 第 ${run} 次运行失败"
#       exit_status=1
#     fi
#   done
# done

# 执行单个任务,多节点
if process_percentile $path ; then
  echo "任务成功: p${p}_${num}"
else
  echo "任务失败: p${p}_${num}" >&2
  exit 1
fi

# 最终状态检查
# if [ "${exit_status}" -eq 0 ]; then
#   echo "所有百分位数处理成功完成"
# else
#   echo "部分处理失败，请检查日志："
#   find "${BASE_DIR}" -name "error.log" -exec grep -l "ERROR" {} \;
# fi

# exit ${exit_status}