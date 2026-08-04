#!/bin/bash

module load miniconda/4.9.2 
source activate qiime2-amplicon-2024.2

study_id=""
file_path=""

# ASVs_counts.tsv to ASVs_counts.biom
biom convert \
  -i ${file_path}/results/ASVs_counts.tsv \
  -o ${file_path}/results/ASVs_counts.biom \
  --table-type="OTU table" \
  --to-hdf5

qiime tools import \
  --input-path ${file_path}/results/ASVs.fa \
  --output-path ${file_path}/results/rep-seqs.qza \
  --type 'FeatureData[Sequence]'

qiime tools import \
  --input-path ${file_path}/results/ASVs_counts.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path ${file_path}/results/table.qza
  
qiime vsearch cluster-features-closed-reference \
    --i-sequences ${file_path}/results/rep-seqs.qza \
    --i-table ${file_path}/results/table.qza \
    --i-reference-sequences SILVA_138.2_SSURef_NR99.qza \
    --p-perc-identity 0.97 --p-threads 12 \
    --p-strand both \
    --o-clustered-table ${file_path}/results/clustered-table \
    --o-clustered-sequences ${file_path}/results/clustered-sequences \
    --o-unmatched-sequences ${file_path}/results/unmatched-sequences

qiime tools export \
  --input-path ${file_path}/results/clustered-table.qza \
  --output-path ${file_path}/results/extracted-feature-table
  
