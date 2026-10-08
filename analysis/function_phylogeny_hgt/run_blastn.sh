#!/bin/bash
#SBATCH --job-name=blastn
#SBATCH -N 1
#SBATCH -n 28
#SBATCH --mem=64G
#SBATCH -o blastn.out
#SBATCH -e blastn.err

# Candidate horizontal gene transfer between genome pairs (Methods: BLASTn
# -perc_identity 99 -word_size 28; an alignment of >= 500 bp at >= 99% identity
# counts as one candidate HGT region). Run from analysis/function_phylogeny_hgt:
#
#   GENOME_DIR=/path/to/genomes sbatch run_blastn.sh      # or: bash run_blastn.sh
#
# GENOME_DIR holds the 1,112 representative genomes as nucleotide FASTA, one
# file per genome named <id>.fna, where <id> is a line of data/genome_id_file
# (e.g. GCA_000154345.1_genomic.fna). They are public NCBI assemblies; the
# accession is the part of <id> before "_genomic".
#
# Writes data/blastn_results/filtered/<id1>.fna_vs_<id2>.fna_filtered.tsv for
# every pair in data/genome_pairs_vsearch.txt -- the files get_HGT_result.py
# reads. Needs blastn/makeblastdb (BLAST+ 2.15.0 was used) and GNU parallel on
# PATH; set CONDA_ENV to activate an environment that provides them.
set -euo pipefail

GENOME_DIR="${GENOME_DIR:?set GENOME_DIR to the directory holding <genome_id>.fna}"
THREADS="${SLURM_NTASKS:-${THREADS:-28}}"
DB_DIR=data/databases_blastn
OUT_DIR=data/blastn_results

if [[ -n "${CONDA_ENV:-}" ]]; then
  set +u; eval "$(conda shell.bash hook)"; conda activate "$CONDA_ENV"; set -u
fi
for tool in makeblastdb blastn parallel; do
  command -v "$tool" >/dev/null || { echo "$tool not found on PATH" >&2; exit 1; }
done

mkdir -p "$DB_DIR" "$OUT_DIR/filtered"

# One nucleotide BLAST database per genome.
while read -r id; do
  [[ -e "$DB_DIR/$id.nsq" || -e "$DB_DIR/$id.00.nsq" ]] && continue
  makeblastdb -in "$GENOME_DIR/$id.fna" -dbtype nucl -out "$DB_DIR/$id" > /dev/null
done < data/genome_id_file

# One line of the pair list ("<id1> <id2>"): query id1 against id2's database.
process_pair() {
  local q s raw
  read -r q s <<< "$1"
  raw="$OUT_DIR/$q.fna_vs_$s.fna.tsv"
  blastn -query "$GENOME_DIR/$q.fna" -db "$DB_DIR/$s" \
         -outfmt "6 qseqid sseqid pident length qstart qend sstart send" \
         -perc_identity 99 -word_size 28 -num_threads 1 -out "$raw"
  awk -F'\t' '$4 >= 500 && $3 >= 99' "$raw" > "$OUT_DIR/filtered/$q.fna_vs_$s.fna_filtered.tsv"
}
export -f process_pair
export GENOME_DIR DB_DIR OUT_DIR

parallel -j "$THREADS" process_pair :::: data/genome_pairs_vsearch.txt
