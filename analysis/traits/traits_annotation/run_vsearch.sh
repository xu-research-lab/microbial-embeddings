#!/bin/bash

module load miniconda/4.9.2 
source activate qiime2-amplicon-2024.2

queries="../data/feces_seq_16S_silva.fasta"
database="../data/99_otus.fasta"
vsearch --usearch_global ${queries} --db ${database} --id 0.99 --blast6out ../data/vsearch.out
