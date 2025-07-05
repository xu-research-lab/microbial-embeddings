#!/bin/python
import sys
import os
import pandas as pd
from carveme import config, project_dir
from reframed.io.sbml import sanitize_id
from carveme.reconstruction.diamond import run_blast, load_diamond_results


def build_model_id(name):
    model_id = sanitize_id(name)
    if not model_id[0].isalpha():
        model_id = 'm_' + model_id
    return model_id

genome_id = sys.argv[1]

inputfile = f"Data/genome/{genome_id}.faa"
model_id = os.path.splitext(os.path.basename(inputfile))[0]
model_id = build_model_id(model_id)
diamond_db = project_dir + config.get('generated', 'diamond_db')
blast_output = f"Data/blast_output/{genome_id}.tsv"
exit_code = run_blast(inputfile, 'protein', blast_output, diamond_db)
annotations = load_diamond_results(blast_output)
gene2gene = annotations.sort_values(by='score', ascending=False) \
                          .groupby('BiGG_gene', as_index=False).apply(lambda x: x.iloc[0])
gene2gene.to_csv(blast_output, index=None)


