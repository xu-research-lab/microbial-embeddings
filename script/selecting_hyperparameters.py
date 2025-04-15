import sys
import biom
import random
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_array

glove_NLP = "/home/dongbiao/word_embedding_microbiome/programe_test/x-lang/GloVe-master/cooccurrence.bin"
data_dir = "/home/dongbiao/word_embedding_microbiome/programe_test/AGP_run_20220702"
cooccurrence_name = ["GloVe", "AGP_russell_rao", "AGP_russell_rao_weight", 
                     "AGP_jaccard", "AGP_faith", "AGP_abundance-percentile",
                     "AGP_abundance-totalsum", "AGP_braycurtis_percentile",
                     "AGP_braycurtis_normalize"]
X_max = [100, 0.1, 0.01, 0.2, 0.04, 0.05, 0.001, 0.1, 0.1]

def plot_xmax_distribution(cooccurrence):
    dt = np.dtype([('word1', 'int32'), ('word2', 'int32'),
                   ('value', 'float64')])
    table = np.fromfile(cooccurrence, dtype=dt)
    
    table = [i[2] for i in table]
    
    p = sum(np.array(table) < 100) / len(table)
    
    return table, p

fig, ax = plt.subplots(3, 3, figsize=(20,20))
n = 0 
for i in range(0, 3):
    for j in range(0, 3):
        if n == 0:
            cooccurrence = glove_NLP
        else:
            cooccurrence = data_dir + "/" + cooccurrence_name[n] + ".bin"
        table, p = plot_xmax_distribution(cooccurrence)
        hist, edges = np.histogram(random.sample(table, 10000), bins=50)
        freq = hist / float(hist.sum())
        width = np.diff(edges)
        ax[i][j].bar(edges[1:], np.log10(freq), width=width, align="edge", ec="k")
        ax[i][j].axvline(x=X_max[n], color='red')
        ax[i][j].set_ylabel("log10(Frequency)")
        ax[i][j].set_title('Co-occurrence distribution of {cooccurrence_name[n]} \n p(co-occurrence < 100)=0.998')
        
        n += 1
fig.savefig('/home/dongbiao/word_embedding_microbiome/programe_test/glove/select_x_max.png', dpi=300)
