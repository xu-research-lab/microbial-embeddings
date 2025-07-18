#!/bin/bash
input="Data/pick_otu.fasta"
aligned="Data/aligned.fasta"
identity_matrix="Data/identity_matrix.txt"
distance_matrix="Data/distance_matrix.txt"

mafft --auto --thread 18 ${input} > ${aligned}


clustalo \
  -i ${aligned} \
  --percent-id \
  --distmat-out=${identity_matrix} \
  --full \
  --force --threads 18 -t DNA