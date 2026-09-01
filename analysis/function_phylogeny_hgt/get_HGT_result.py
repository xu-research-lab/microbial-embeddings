import pandas as pd
import numpy as np

vsearch_res = pd.read_csv("../resources/genome_mapping/data/vsearch_res.csv")
genome_16S_fid = {}
for i in range(vsearch_res.shape[0]):
    temp = vsearch_res.iloc[i].values
    genome_16S_fid[temp[14]] = temp[0]

with open("data/identity_matrix.txt", "r") as file:
    data = file.readlines()
value = []
genome_id = []
for line in data:
    genome_id.append(line.split(" ")[0])
    temp = line.split("\n")[0].split(" ")[1:]
    temp = [float(item) for item in temp if item != ""]
    value.append(temp)
genome_id = [genome_16S_fid.get(i, i) for i in genome_id]
identity_table = pd.DataFrame(data=value[1:], index=genome_id[1:], columns=genome_id[1:])
identity_table = identity_table.loc[vsearch_res.query_id.values, vsearch_res.query_id.values]

genome_fid = {}
for i in range(vsearch_res.shape[0]):
    temp = vsearch_res.iloc[i].values
    genome_fid[temp[15]] = temp[0]

genome_pairs = pd.read_csv("data/genome_pairs_vsearch_new.txt", header=None,  sep=" ")

path="data/blastn_results/filtered"
hgt = []
genome_1 = genome_pairs.iloc[:,0].values
genome_2 = genome_pairs.iloc[:,1].values
for i in tqdm(range(genome_pairs.shape[0]), desc="Processing"):
    temp = genome_pairs.iloc[i].values
    try:
        table = pd.read_csv(f"{path}/{temp[0]}.fna_vs_{temp[1]}.fna_filtered.tsv")
        hgt.append(table.shape[0])
    except Exception:
        hgt.append(0)

identity = []
cosine_co = []
phy_dis = []
for i in range(len(genome_1)):
    id_1 = genome_fid[genome_1[i]]
    id_2 = genome_fid[genome_2[i]]
    identity.append(identity_table.loc[id_1, id_2])
    cosine_co.append(embed_cos.loc[id_1, id_2])
    phy_dis.append(dm.loc[id_1, id_2])

fid_1 = [genome_fid[i] for i in genome_1]
fid_2 = [genome_fid[i] for i in genome_2]

hgt_res = pd.DataFrame({"id_1":fid_1, "id_2":fid_2, "identity":identity, 
                        "cosine_co":cosine_co, "hgt":hgt, "phy_dis":phy_dis})

hgt_res.to_csv("/home/dongbiao/word_embedding_microbiome/HGT/results/hgt.csv", index=None)

def min_max_normalization(data): 
    data_min = np.min(data)
    data_max = np.max(data)
    normalized_data = (data - data_min) / (data_max - data_min)
    return normalized_data

with open("/home/dongbiao/word_embedding_microbiome/HGT/results/genome_fid.json", "r") as f:
    genome_fid = json.load(f) 
phy_dis_vector = np.linspace(0.02, 1, 10)
embed_dis_vector = np.linspace(0.1, 0.80, 10)
