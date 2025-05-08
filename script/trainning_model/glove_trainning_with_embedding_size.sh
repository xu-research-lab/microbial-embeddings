#!/bin/bash

#SBATCH --job-name=glove_embeding_size
#SBATCH -N 1  
#SBATCH -p cu  
#SBATCH -n 28                 # 使用28核
#SBATCH --mem=250G  
#SBATCH -o master_%a_r.log  # 改动点
#SBATCH -e master_%a_r.err  # 改动点
#SBATCH --array=1-4%2
#SBATCH --exclude=cu01

# 激活conda环境
source /home/cjj/miniconda3/etc/profile.d/conda.sh
conda activate microbiome_deep || {
  echo "无法激活conda环境"
  exit 1
}
which python


# 读取参数矩阵
mapfile -t params < matrix_newdata_ab_embedsize.txt  # TODO 改动点
IFS=' ' read -r p num embedding_size <<< "${params[$SLURM_ARRAY_TASK_ID-1]}"

# id=$(sed "${SLURM_ARRAY_TASK_ID}q;d" matrix.txt)

# 全局配置
GENERATE_COOC=0                  # 0=复用现有文件，1=生成新文件 # 改动点
SOURCE_DIR="p80_abundance_percentile"  # 复用模式时源目录
BASE_DIR="."           # 输出基准目录
INPUT_BIOM="${BASE_DIR}/table_permutation.biom"  # TODO 输入文件路径
CPUS_PER_TASK=28                       # 每个任务CPU数



# 处理单个百分位数的函数
process_percentile() {
  # TODO 固定距离调参 80
  local percentile=$1
  local run=$2
  local embedding_size="$3"
  local work_dir="${BASE_DIR}/p${percentile}_${run}_emb${embedding_size}" # 改动点
  local result_dir="${work_dir}/result"
  local log_dir="${work_dir}/log"

  # TODO 距离调参
  local metric=${run}

  rm -rf "${result_dir}" "${log_dir}"
  mkdir -p "${work_dir}" "${result_dir}" "${log_dir}"

  echo "==============================================="
  echo "[$(date)] 开始处理百分位数 ${percentile}%"
  echo "工作目录: ${work_dir}"
  echo "日志目录: ${log_dir}"
  
  # 创建工作目录并复制biom文件
  if [[ ! -f "${INPUT_BIOM}" ]]; then
    echo "错误：输入文件不存在 ${INPUT_BIOM}" | tee -a "${log_dir}/error.log"
    return 1
  fi
  cp "${INPUT_BIOM}" "${work_dir}/input.biom" || return 1

  # 0=复用现有文件，1=重新生成新文件共现向量
  if [[ "${GENERATE_COOC}" -eq 1 ]]; then

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
  else
    # ******** 复用模式 ********
    echo "复用现有文件从: ${SOURCE_DIR}" | tee -a "${log_dir}/process.log"
    local required_files=(
      "${SOURCE_DIR}/feature-dict.csv"
      "${SOURCE_DIR}/table.co" 
    )
    
    # 检查源文件是否存在
    for f in "${required_files[@]}"; do
      if [[ ! -f "${f}" ]]; then
        echo "缺失必要文件: ${f}" | tee -a "${log_dir}/error.log"
        return 1
      fi
      cp "${f}" "${work_dir}/" || return 1
    done
  fi

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
    --embedding-size ${embedding_size} \
    --iter 100 \
    --cpus ${CPUS_PER_TASK} \
    > "${log_dir}/step3.log" 2>&1 || {
    echo "模型训练失败！百分位数: ${percentile}" | tee -a "${log_dir}/error.log"
    return 1
  }

  echo "[$(date)] 百分位数 ${percentile}% 处理完成" | tee -a "${log_dir}/process.log"
  return 0
}

# 执行单个任务,多节点
if process_percentile $p $num $embedding_size; then
  echo "任务成功: p${p}_${num}_${embedding_size}"
else
  echo "任务失败: p${p}_${num}_${embedding_size}" >&2
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