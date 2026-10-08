#!/bin/bash
#SBATCH --job-name=mafft
#SBATCH -N 1
#SBATCH -n 18
#SBATCH --mem=50G
#SBATCH -o mafft.out
#SBATCH -e mafft.err

# Pairwise 16S identity of the representative genomes (Methods: MAFFT --auto,
# then Clustal Omega --percent-id --full). Run from analysis/function_phylogeny_hgt.
# Needs mafft and clustalo on PATH; set CONDA_ENV to activate an environment that
# provides them (e.g. conda create -n hgt_align -c bioconda mafft=7.525 clustalo=1.2.4).
set -euo pipefail

THREADS="${SLURM_NTASKS:-${THREADS:-18}}"
if [[ -n "${CONDA_ENV:-}" ]]; then
  set +u; eval "$(conda shell.bash hook)"; conda activate "$CONDA_ENV"; set -u
fi
for tool in mafft clustalo; do
  command -v "$tool" >/dev/null || { echo "$tool not found on PATH" >&2; exit 1; }
done

input="data/pick_otu.fasta"
aligned="data/aligned.fasta"
identity_matrix="data/identity_matrix.txt"

mafft --auto --thread "$THREADS" "$input" > "$aligned"

clustalo \
  -i "$aligned" \
  --percent-id \
  --distmat-out="$identity_matrix" \
  --full \
  --force --threads "$THREADS" -t DNA
