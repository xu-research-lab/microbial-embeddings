#!/bin/bash

study_id=""
file_path=""
path=""

# makedirs
mkdir ${file_path}/intermediate ${file_path}/rawdata ${file_path}/results ${file_path}/temp

# download rawdata
prefetch --option-file ${file_path}/sample.txt --output-directory ${file_path}/sra
ls ${file_path}/sra/*/*.sra | awk -F'/' '{print $11}' | awk -F'.' '{print $1}' > ${file_path}/sample.txt

threads=20
# fastq-dump
cat "${file_path}/sample.txt" | ${path}/parallel -j ${threads} --bar \
    "fastq-dump --split-3 ${file_path}/sra/{}/{}.sra --outdir ${file_path}/rawdata"

# Run DADA2
ls ${file_path}/rawdata/ > ${file_path}/run_sample.txt

Rscript process_forwards.R ${file_path} trim
