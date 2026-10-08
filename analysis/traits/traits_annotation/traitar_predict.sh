#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem 250G

# Traitar phenotype prediction for the representative genomes (Methods: Traitar
# v1.1.1, default parameters). Run from analysis/traits/traits_annotation:
#
#   PFAM_DIR=/path/to/pfam GENOME_DIR=/path/to/genomes sbatch traitar_predict.sh
#
# PFAM_DIR   : the Pfam database Traitar uses (download once with `traitar pfam`;
#              Pfam 33.1 was used)
# GENOME_DIR : the genome nucleotide FASTA files listed in the sample_file_name
#              column of samples.txt (e.g. GCA_000154345.1_genomic.fna)
# OUT_DIR    : output directory (default traitar_out); its phenotype table is
#              what traits_predict_Traitar.csv was built from.
# Needs `traitar` on PATH; set CONDA_ENV to activate an environment that has it.
set -euo pipefail

PFAM_DIR="${PFAM_DIR:?set PFAM_DIR to the Pfam directory Traitar uses}"
GENOME_DIR="${GENOME_DIR:?set GENOME_DIR to the directory holding the genome FASTA files}"
OUT_DIR="${OUT_DIR:-traitar_out}"
CPUS="${SLURM_CPUS_PER_TASK:-28}"

if [[ -n "${CONDA_ENV:-}" ]]; then
  set +u; eval "$(conda shell.bash hook)"; conda activate "$CONDA_ENV"; set -u
fi
command -v traitar >/dev/null || { echo "traitar not found on PATH" >&2; exit 1; }

traitar phenotype "$PFAM_DIR" "$GENOME_DIR" samples.txt from_nucleotides "$OUT_DIR" \
    -o -c "$CPUS" -x 14
