import pandas as pd
import random
import numpy as np
import biom
import pickle
from sklearn.metrics.pairwise import cosine_similarity
from skbio.stats.distance import permanova
from skbio import DistanceMatrix
from skbio.stats.distance import anosim

AGP_table = biom.load_table("/home/dongbiao/word_embedding_microbiome/Classification_prediction/data/AGP_Glove.biom")
sid_1 = AGP_table.ids(axis="observation")
all_feces_table = biom.load_table("/home/dongbiao/word_embedding_microbiome/databaes/16Sdatabase/human16S_feces/human_train.biom")
sid_2 = all_feces_table.ids(axis="observation")
sid = np.intersect1d(sid_1, sid_2)
sid_sample = random.sample(sid.tolist(), 5000)

# 代谢通路表
pathway_table = pd.read_csv("/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/pathway/OTU_pathway_table.csv",
                            index_col=0)
                            
sid_sample = np.intersect1d(sid_sample, pathway_table.index.tolist())
pathway_table = pathway_table.loc[[int(i) for i in sid_sample], ]
total_col = pathway_table.shape[0]
keep1 = np.array(pathway_table.sum()) < round(total_col * 0.9, 0)
keep2 = np.array(pathway_table.sum()) > round(total_col * 0.1, 0)
pathway_table = pathway_table.loc[:, keep1 & keep2]

otu_embedding = ["russell_rao", "russell_rao_weight", "faith", "jaccard", 
                 "abundance-percentile", "abundance-totalsum", 
                 "braycurtis-totalsum", "braycurtis-percentile", "phylogeny", "PCA"]

embedding_dir = "/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702"
embedding_type = "agp"

res = {}
for j in otu_embedding:
    embedding_vector = pd.read_csv(f"{embedding_dir}/{j}/{j}_100.txt",
                                   sep=" ", index_col=0, header=None, dtype={0:str})
    embedding_vector = embedding_vector.loc[sid_sample, ]
    distance = cosine_similarity(embedding_vector)
    
    for i in range(len(distance)):
        distance[i, i] = 0
    
    dm = DistanceMatrix(distance)
    R_statistic = []
    for i in range(pathway_table.shape[1]):
        grouping = pathway_table.iloc[:, i].tolist()
        tmp = anosim(dm, grouping, permutations=99)
        R_statistic.append(tmp['test statistic'])

    res[f'{j}'] = R_statistic

otuput_dir = "/home/dongbiao/word_embedding_microbiome/programe_test/glove/glove_embedding/"
with open(f'{otuput_dir}/R_res_{embedding_type}.pkl', 'wb') as f:
    pickle.dump(res, f)

