#!/bin/python
import pandas as pd
import numpy as np
from itertools import combinations
from reframed import load_cbmodel

from joblib import Parallel, delayed
from tqdm import tqdm
from reframed import set_default_solver

def jaccard_metabolic_dism(model_1, model_2):
    
    model1 = load_cbmodel(f"/home/dongbiao/word_embedding_microbiome/modelseed/predict_otu_metabolic/OTU_metabolic_model_WD/{model_1}.xml")
    model2 = load_cbmodel(f"/home/dongbiao/word_embedding_microbiome/modelseed/predict_otu_metabolic/OTU_metabolic_model_WD/{model_2}.xml")
    
    # 获取两个模型中的反应集合
    reactions_model1 = set([r for r in model1.reactions])
    reactions_model2 = set([r for r in model2.reactions])
    
    # 计算 Jaccard 距离
    intersection = reactions_model1.intersection(reactions_model2)  # 交集
    union = reactions_model1.union(reactions_model2)  # 并集
    jaccard_distance = 1 - len(intersection) / len(union)

    return jaccard_distance


# 读取数据
run_metabolic_data = pd.read_csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/run_metabolic_data_module_33.csv")
# 定义函数来处理每一行数据
def process_row(row):
    temp = row.values
    jaccard_dist = jaccard_metabolic_dism(temp[0], temp[1])
    community = temp[2]
    type_ = temp[3]
    return jaccard_dist, community, type_

# 使用 joblib 并行处理
results = Parallel(n_jobs=20)(delayed(process_row)(run_metabolic_data.iloc[i]) for i in tqdm(range(run_metabolic_data.shape[0]), desc="Processing"))

# 将结果拆分为各自的列表
jaccard_distance, community, types = zip(*results)

# 如果需要将结果转换为列表
jaccard_distance = list(jaccard_distance)
community = list(community)
types = list(types)

metabolic_dis_res = pd.DataFrame({"community":community, "jaccard_distance":jaccard_distance, "group":types})
metabolic_dis_res.to_csv("/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/metabolic_dissimilarity_module_33.csv", index=None)