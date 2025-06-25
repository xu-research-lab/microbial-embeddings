
from cobra.io import read_sbml_model
from tqdm import tqdm
from joblib import Parallel, delayed
import pandas as pd
import pickle  

co_clusters = pd.read_csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/cluster/spectral_clustering.csv")
co_clusters.fid = [str(i) for i in co_clusters.fid]

all_fid = co_clusters.fid.values
# set_default_solver("cplex")
# 计算最小培养基
# 默认目标是最小化培养基中的营养物质总量
# 定义一个函数来处理每个 i 的任务
def process_model(i):
    # 加载模型
    model = read_sbml_model(f"/home/dongbiao/word_embedding_microbiome/modelseed/predict_otu_metabolic/OTU_metabolic_model_WD/{i}.xml")
    try:
        model.medium.pop("EX_o2_e")
    except Exception as e:
        print("No EX_o2_e")
    solution = model.optimize()
    secretions = []
    for rxn in model.exchanges:
        flux = solution.fluxes[rxn.id]
        if flux > 1e-6:  # 忽略微小通量
            metabolite = rxn.reactants[0]  # 获取代谢物对象
            secretions.append(metabolite.name)

    return i, secretions

# 使用 joblib 并行处理
results = Parallel(n_jobs=10)(delayed(process_model)(i) for i in tqdm(all_fid, desc="Processing"))

# 将结果转换为字典
compounds_comptete = {i: compounds for i, compounds in results}

with open('/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/compounds_compete_produce.pkl', 'wb') as file:  
    pickle.dump(compounds_comptete, file) 
