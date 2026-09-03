#!/bin/bash
source activate picrust2
hsp.py -t data/bac.tre --observed_trait_table mapping_bigg_gene_table.tsv -o bigg_gene_predicted.tsv -p 4 -m mp -n
hsp.py -t data/bac.tre --observed_trait_table mapping_scores_table.tsv -o scores_predicted.tsv -p 24 -m pic -n