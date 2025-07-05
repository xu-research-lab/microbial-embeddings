#!/bin/python

import os
import subprocess
from smetana.interface import main
import pandas as pd
import numpy as np
import sys
from reframed import set_default_solver
from joblib import Parallel, delayed
from tqdm import tqdm

mediadb = "Data/media_db.tsv"

def maincall(i, run_file):
    
    os.mkdir(f"{run_file}/{i}")
    run_file = f"{run_file}/{i}"
    genome_id = modelseed_micro_pairs.loc[i, ].values
    split_id = i
    metabolic_model_path = "Data/OTU_metabolic_model_M3"
    output = "Data/smetana/results"
    
    ### generate community file
    id_2 = [f"{i}_{split_id}" for i in genome_id]
    id_1 = ["comm"] * len(id_2)
    coum = pd.DataFrame({"id_1":id_1, "id_2":id_2})
    coum.to_csv(f"{run_file}/communities_{split_id}.tsv", header=None, index=None, sep="\t")
    for j in range(len(genome_id)):
        subprocess.run(['cp', f'{metabolic_model_path}/{genome_id[j]}.xml', f'{run_file}/{genome_id[j]}_{split_id}.xml'], 
                        capture_output=True, text=True)
    otuput_file = '_'.join(genome_id)
    
    ### model: global, detailed
    main([f"{run_file}/*_{split_id}.xml"], mode="global", output=f"{output}/{otuput_file}_WD_output", media="M11",
        mediadb=mediadb, exclude="Data/inorganic.txt", communities=f"{run_file}/communities_{split_id}.tsv", 
        use_lp=True, ignore_coupling=True)
    

    
    subprocess.run(['rm', '-rf', f'{run_file}'], capture_output=True, text=True)
    
if __name__ == '__main__':

    set_default_solver("cplex")
    input_file = sys.argv[1]
    run_file = sys.argv[2]
    modelseed_micro_pairs = pd.read_csv(input_file, dtype=str)
    Parallel(n_jobs=6)(delayed(maincall)(i, run_file) for i in tqdm(range(modelseed_micro_pairs.shape[0]), desc="Processing"))
