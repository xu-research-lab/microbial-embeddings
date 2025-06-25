
from reframed import load_cbmodel, minimal_medium
from reframed.io.sbml import load_cbmodel
from joblib import Parallel, delayed
from reframed import set_default_solver
from tqdm import tqdm
import pandas as pd
import pickle  

co_clusters = pd.read_csv("/home/dongbiao/word_embedding_microbiome/programe_test/embedding_explain/cluster/spectral_clustering.csv")
co_clusters.fid = [str(i) for i in co_clusters.fid]

# module_93 = co_clusters.loc[co_clusters.clusters == 93].fid.unique()
module_111 = co_clusters.loc[co_clusters.clusters == 111].fid.unique()
# all_fid = co_clusters.fid.values
set_default_solver("cplex")
# 计算最小培养基
# 默认目标是最小化培养基中的营养物质总量
# 定义一个函数来处理每个 i 的任务
def process_model(i):
    # 加载模型
    model = load_cbmodel(f"/home/dongbiao/word_embedding_microbiome/modelseed/predict_otu_metabolic/OTU_metabolic_model_WD/{i}.xml")
    # 移除反应
    model.remove_reaction('R_EX_o2_e')
    # 计算最小培养基
    min_medium = minimal_medium(model, min_growth=0.1)  # min_growth 是最小生长速率
    # 返回结果
    return i, list(min_medium[0])

# 使用 joblib 并行处理
results = Parallel(n_jobs=1)(delayed(process_model)(i) for i in tqdm(module_111, desc="Processing"))

# 将结果转换为字典
compounds_comptete = {i: compounds for i, compounds in results}

 
with open('/home/dongbiao/word_embedding_microbiome/modelseed/metabolic_dissimilarity/compounds_compete_module_111.pkl', 'wb') as file:  
    pickle.dump(compounds_comptete, file) 
