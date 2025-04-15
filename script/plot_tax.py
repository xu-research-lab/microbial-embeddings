import pandas as pd
import random
from sklearn.manifold import TSNE
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import biom


taxonomy = pd.read_table("/beegfs/db/greengenes/gg_13_8_otus/taxonomy/99_otu_taxonomy.txt", header=None, sep='\t')
tax = taxonomy[1].str.split(';',expand=True).add_prefix('taxonomy_')
tax_id = list(taxonomy[0])

AGP_table = biom.load_table("/home/dongbiao/word_embedding_microbiome/Classification_prediction/data/AGP_Glove.biom")
sid_1 = AGP_table.ids(axis="observation")
all_feces_table = biom.load_table("/home/dongbiao/word_embedding_microbiome/databaes/16Sdatabase/human16S_feces/human_train.biom")
sid_2 = all_feces_table.ids(axis="observation")
sid = np.intersect1d(sid_1, sid_2)
sid_sample = random.sample(sid.tolist(), 5000)

def plot_tsne(embedding_vector, taxonomy, num, level, title, save_path, co_method, type_embedding):
    distance = cosine_similarity(embedding_vector)
    tsne = TSNE(n_components=2)
    result = tsne.fit_transform(distance)
    
    fid_ = [int(i) for i in sid_sample]
    id = [tax_id.index(i) for i in fid_]
    tmp = tax.iloc[id].groupby([f'{taxonomy}']).count().sort_values(by="taxonomy_0", ascending=False)[0:num].reset_index()
    
    df = pd.DataFrame(dict(x=result[:, 0], y=result[:, 1], color=tax.iloc[id][f'{taxonomy}']))
    df = pd.merge(df, tmp, how='left', left_on = 'color', right_on = f'{taxonomy}').replace(np.nan, 'others')

    color_labels = df[f'{taxonomy}'].unique()

    # List of RGB triplets
    rgb_values = sns.color_palette("tab20", num+1)

    # Finally use the mapped values
    plt.figure(figsize=(7,5))
    g = sns.scatterplot(x='x', y = 'y', s=30, data=df, hue = f'{taxonomy}', palette=rgb_values, sizes=100)
    box = g.get_position()
    g.set_position([box.x0, box.y0, box.width * 0.8, box.height]) 
    plt.legend(title = f'{level}', fontsize = 8, title_fontsize = 8, 
               loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title(f"tSNE of {title} embedding of OTUs")
    plt.grid()
    plt.savefig(f"{save_path}/tsne_{co_method}_{level}_{type_embedding}.png", dpi=100)
    plt.close()
    
    
otu_embedding = ["russell_rao", "russell_rao_weight", "faith", "jaccard", 
                 "abundance-percentile", "abundance-totalsum", 
                 "braycurtis-totalsum", "braycurtis-percentile", "phylogeny", "PCA"]


n=0
for i in otu_embedding:
    names = ["Russell_rao", "Russell_rao_weight", "Faith", "Jaccard", 
             "Abundance_percentile", "Abundance_totalsum", 
             "Braycurtis_totalsum", "Braycurtis_percentile", "phylogeny_glove", "PCA"]
    qual_vec_file_1 = f"/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702/{i}/{i}_100.txt"
    qual_vec_file_2 = f"/home/dongbiao/word_embedding_microbiome/databaes/16Sdatabase/glove_1/{i}/{i}_100.txt"

    # AGP
    qual_vecs = pd.read_csv(qual_vec_file_1, sep=" ", index_col=0, header=None, dtype={0:str})
    qual_vecs = qual_vecs.loc[sid_sample, ]
    plot_tsne(qual_vecs, 'taxonomy_1', 17, 'Phylum', names[n], 
              "/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/res/tax_plot", i, "agp")
    plot_tsne(qual_vecs, 'taxonomy_5', 17, 'Genes', names[n], 
              "/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/res/tax_plot", i, "agp")

    # all feces
    # qual_vecs = pd.read_csv(qual_vec_file_2, sep=" ", index_col=0, header=None, dtype={0:str})
    # qual_vecs = qual_vecs.loc[sid_sample, ]
    # plot_tsne(qual_vecs, 'taxonomy_1', 17, 'Phylum', names[n], 
    #           "/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/res/tax_plot", i, "all_feces")
    # plot_tsne(qual_vecs, 'taxonomy_5', 17, 'Genes', names[n], 
    #           "/home/dongbiao/word_embedding_microbiome/programe_test/phylogeny/res/tax_plot", i, "all_feces")
    n += 1



    
