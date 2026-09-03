#!/usr/bin/env bash
#SBATCH --job-name=sne_glove
#SBATCH --nodes=1
#SBATCH --partition=cu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=250G
#SBATCH --output=results/slurm_%j.out
#SBATCH --error=results/slurm_%j.err
#SBATCH --exclude=cu01

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "$script_dir/../../.." && pwd)

biom_file=${1:-$repo_root/data/gut_pretraining.biom}
metric=${2:-abundance_percentile}
embedding_size=${3:-100}
output_dir=${4:-$script_dir/results/${metric}_${embedding_size}}
cpus=${5:-${SLURM_CPUS_PER_TASK:-28}}
cooccurrence_dir=${6:-}
xmax_percentile=${7:-80}

[[ "$biom_file" = /* ]] || biom_file="$repo_root/$biom_file"
[[ "$output_dir" = /* ]] || output_dir="$repo_root/$output_dir"
if [[ -n "$cooccurrence_dir" && "$cooccurrence_dir" != /* ]]; then
  cooccurrence_dir="$repo_root/$cooccurrence_dir"
fi

if [[ ! -f "$biom_file" ]]; then
  echo "Input BIOM file not found: $biom_file" >&2
  exit 1
fi
if [[ -e "$output_dir" ]]; then
  echo "Output already exists: $output_dir" >&2
  exit 1
fi

source /home/cjj/miniconda3/etc/profile.d/conda.sh
conda activate microbiome_deep
export PYTHONPATH="$repo_root${PYTHONPATH:+:$PYTHONPATH}"
export OMP_NUM_THREADS="$cpus"

mkdir -p "$output_dir"
if [[ -n "$cooccurrence_dir" ]]; then
  feature_dict="$cooccurrence_dir/feature-dict.txt"
  [[ -f "$feature_dict" ]] || feature_dict="$cooccurrence_dir/feature-dict.csv"
  cp "$feature_dict" "$output_dir/feature-dict.txt"
  cp "$cooccurrence_dir/table.co" "$output_dir/table.co"
else
  membed dict -b "$biom_file" -d "$output_dir/feature-dict.txt"
  membed cooccur -b "$biom_file" -c "$output_dir/table.co" --metric "$metric" --cpus "$cpus"
fi
membed build-x-max-file -c "$output_dir/table.co" -x "$output_dir/xmax.npy" --percentile_num "$xmax_percentile"
membed glove-train \
  -d "$output_dir/feature-dict.txt" \
  -c "$output_dir/table.co" \
  -x "$output_dir/xmax.npy" \
  -r "$output_dir" \
  --lr 0.05 \
  --embedding-size "$embedding_size" \
  --iter 100 \
  --cpus "$cpus"
