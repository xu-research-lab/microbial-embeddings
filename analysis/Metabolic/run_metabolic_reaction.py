#!/bin/python
import pandas as pd
import numpy as np
from itertools import combinations
from reframed import load_cbmodel
import pickle  
from joblib import Parallel, delayed
from tqdm import tqdm

co_clusters = pd.read_csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/cluster/spectral_clustering.csv")
co_clusters.fid = [str(i) for i in co_clusters.fid]
all_fid = co_clusters.fid.values


def process_model(i):
    # 加载模型
    model = load_cbmodel(f"/home/dongbiao/word_embedding_microbiome/modelseed/predict_otu_metabolic/OTU_metabolic_model_WD/{i}.xml")
    reactions_model = [r for r in model.reactions]
    
    return i, reactions_model

# 使用 joblib 并行处理
results = Parallel(n_jobs=10)(delayed(process_model)(i) for i in tqdm(all_fid, desc="Processing"))

# 将结果转换为字典
compounds_comptete = {i: reactions_model for i, reactions_model in results}

with open('/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/metabolic_reaction.pkl', 'wb') as file:  
    pickle.dump(compounds_comptete, file) 
