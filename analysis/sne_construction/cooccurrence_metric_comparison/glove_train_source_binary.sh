#!/bin/bash

#SBATCH --job-name=glove_array
#SBATCH -N 1
#SBATCH -p cu
#SBATCH -n 28                  # 使用28核
#SBATCH --mem=100G
#SBATCH -o log/%x_%a_b.out       # Slurm 标准输出日志
#SBATCH -e log/%x_%a_b.err       # Slurm 错误输出日志
#SBATCH --array=1-21%3           # 提交40个任务 (5个数据集 * 8个度量)
#SBATCH --exclude=cu01

# --- 全局配置 ---
BASE_OUTPUT_DIR="./" # 所有结果的总输出目录
PARAM_FILE="param_matrix_binary.txt"    # 参数矩阵文件
CPUS_PER_TASK=28                 # 每个任务使用的CPU核心数
PERCENTILE_FOR_XMAX=80           # 用于 'build-x-max-file' 的百分位数值

# --- 环境激活 ---
echo "加载 Conda 环境..."
source /home/cjj/miniconda3/etc/profile.d/conda.sh
conda activate microbiome_deep || {
  echo "错误：无法激活 Conda 环境 'microbiome_deep'"
  exit 1
}
echo "Python anaconda环境: $(which python)"
echo "==============================================="


# --- 读取当前任务的参数 ---
# 从参数矩阵文件中读取第 $SLURM_ARRAY_TASK_ID 行
mapfile -t params < ${PARAM_FILE}
IFS=' ' read -r biom_path metric <<< "${params[$SLURM_ARRAY_TASK_ID-1]}"

# 从biom路径中提取数据集名称，用于创建目录
dataset_name=$(basename "$biom_path" .biom)


# --- 设置工作目录和日志目录 ---
# 根据度量和数据集名称创建唯一的输出目录
work_dir="${BASE_OUTPUT_DIR}/${metric}/${dataset_name}"
result_dir="${work_dir}/result_vectors" # 存放最终的嵌入向量
job_log_dir="${work_dir}/logs"          # 存放此任务的详细步骤日志
local_biom_path="${work_dir}/input.biom"  # 工作目录中biom文件的路径

# 清理并创建目录
rm -rf "${work_dir}"
mkdir -p "${result_dir}" "${job_log_dir}"


# --- 开始处理任务 ---
echo "[$(date)] ==> 开始处理任务 #${SLURM_ARRAY_TASK_ID}"
echo "  - 原始数据集: ${biom_path}"
echo "  - 度量      : ${metric}"
echo "  - 工作目录  : ${work_dir}"
echo "==============================================="

# 【新增】拷贝biom文件到工作目录
echo "步骤 0/5: 拷贝 .biom 文件到工作目录..." | tee -a "${job_log_dir}/process.log"
if [[ ! -f "${biom_path}" ]]; then
  echo "错误：输入文件不存在 ${biom_path}" | tee -a "${job_log_dir}/error.log"
  exit 1
fi
cp "${biom_path}" "${local_biom_path}" || {
    echo "错误：拷贝biom文件失败！" | tee -a "${job_log_dir}/error.log"
    exit 1
}

# 步骤 1: 生成特征字典
echo "步骤 1/5: 生成特征字典..." | tee -a "${job_log_dir}/process.log"
membed dict -b "${local_biom_path}" -d "${work_dir}/feature-dict.csv" \
  > "${job_log_dir}/step1_dict.log" 2>&1 || {
  echo "错误：特征字典生成失败！" | tee -a "${job_log_dir}/error.log"
  exit 1
}

# 步骤 2: 生成共现矩阵
echo "步骤 2/5: 生成共现矩阵 (Metric: ${metric})..." | tee -a "${job_log_dir}/process.log"
HDF5_USE_FILE_LOCKING=FALSE membed cooccur \
  -b "${local_biom_path}" \
  -c "${work_dir}/table.co" \
  --metric "${metric}" \
  --cpus "${CPUS_PER_TASK}" \
  > "${job_log_dir}/step2_cooccur.log" 2>&1 || {
  echo "错误：共现矩阵生成失败！" | tee -a "${job_log_dir}/error.log"
  exit 1
}

# 步骤 3: 生成 xmax 值文件
echo "步骤 3/5: 生成 xmax 文件 (Percentile: ${PERCENTILE_FOR_XMAX}%)..." | tee -a "${job_log_dir}/process.log"
membed build-x-max-file -c "${work_dir}/table.co" -x "${work_dir}/xmax_file.npy" \
  --percentile_num ${PERCENTILE_FOR_XMAX} \
  > "${job_log_dir}/step3_xmax.log" 2>&1 || {
  echo "错误：xmax 文件生成失败！" | tee -a "${job_log_dir}/error.log"
  exit 1
}

# 步骤 4: 训练 GloVe 模型
echo "步骤 4/5: 训练 GloVe 模型..." | tee -a "${job_log_dir}/process.log"
export OMP_NUM_THREADS=${CPUS_PER_TASK}
membed glove-train -d "${work_dir}/feature-dict.csv" \
  -c "${work_dir}/table.co" \
  -r "${result_dir}" \
  -x "${work_dir}/xmax_file.npy" \
  --lr 0.05 \
  --embedding-size 100 \
  --iter 100 \
  --cpus "${CPUS_PER_TASK}" \
  > "${job_log_dir}/step4_glove_train.log" 2>&1 || {
  echo "错误：模型训练失败！" | tee -a "${job_log_dir}/error.log"
  exit 1
}

# 步骤 5: 完成
echo "步骤 5/5: 清理与完成" | tee -a "${job_log_dir}/process.log"
echo "==============================================="
echo "[$(date)] ==> 任务 #${SLURM_ARRAY_TASK_ID} 处理完成"
echo "结果保存在: ${result_dir}"
echo "==============================================="

exit 0