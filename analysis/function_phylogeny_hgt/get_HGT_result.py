import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.metrics.pairwise import cosine_similarity
import dendropy
from phylodm import PhyloDM

### emb
emb = pd.read_csv("../../data/social_niche_embedding_100.txt", sep=" ", header=None, index_col=0).drop("<unk>")
embed_cos = pd.DataFrame(cosine_similarity(emb), index=emb.index, columns=emb.index)

tree = dendropy.Tree.get_from_path('../../data/SSURefNR99_1200_slv_138_2_subset.tre', schema='newick')
pdm = PhyloDM.load_from_dendropy(tree)
dm = pdm.dm(norm=False)
# phylodm 3.x returns dm() rows/cols in sorted label order, while taxa() is in tree order
labels = sorted(pdm.taxa())
dm = pd.DataFrame(data=dm, index=labels, columns=labels)

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

genome_pairs = pd.read_csv("data/genome_pairs_vsearch.txt", header=None,  sep=" ")

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

hgt_res.to_csv("data/hgt.csv", index=None)

phy_dis_vector = np.linspace(0.00, 1, 10)
embed_dis_vector = np.linspace(0.1, 0.80, 10)

phy_hgt_rate = []
phy_hgt_rate_high_embed = []
phy_hgt_rate_low_embed = []
for i in range(1, 10):
    temp = hgt_res.loc[(hgt_res.phy_dis.values > phy_dis_vector[i-1]) & (hgt_res.phy_dis.values < phy_dis_vector[i])]
    # embed_sim_median = np.median(temp.cosine_co.values)
    embed_sim_median = 0.6
    temp_high_embed = temp.loc[temp.cosine_co.values > embed_sim_median]
    temp_low_embed = temp.loc[temp.cosine_co.values < embed_sim_median]
    temp_high_embed = temp_high_embed.hgt.values
    temp_low_embed = temp_low_embed.hgt.values
    temp = temp.hgt.values
    
    phy_hgt_rate.append(np.sum(temp  > 0) / len(temp) * 100)
    phy_hgt_rate_high_embed.append(np.sum(temp_high_embed  > 0) / len(temp_high_embed ) * 100)
    phy_hgt_rate_low_embed.append(np.sum(temp_low_embed > 0) / len(temp_low_embed) * 100)

embed_hgt_rate = []
for i in range(1, 10):
    temp = hgt_res.loc[(hgt_res.cosine_co.values > embed_dis_vector[i-1]) & (hgt_res.cosine_co.values < embed_dis_vector[i])]
    temp = temp.hgt.values
    embed_hgt_rate.append(np.sum(temp > 0) / len(temp) * 100)


hgt_plot_res = pd.DataFrame({"distance": list(phy_dis_vector[1:]) + list(phy_dis_vector[1:]) + list(phy_dis_vector[1:]) + list(embed_dis_vector[1:]),
                             "hgt_rate": phy_hgt_rate + phy_hgt_rate_high_embed + phy_hgt_rate_low_embed + embed_hgt_rate,
                             "group": ["Phylo"] * len(phy_hgt_rate) + ["Phylo SNEsim > 0.6"] * len(phy_hgt_rate) + ["Phylo SNEsim < 0.6"] * len(phy_hgt_rate) +  ["SNE"] * len(embed_hgt_rate)})

hgt_plot_res.to_csv("data/hgt_plot_res.csv", index=None)