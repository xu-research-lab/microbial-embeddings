import sys
import biom
import dendropy
import numpy as np
from scipy.sparse import coo_array
from phylodm import PhyloDM
import matplotlib.pyplot as plt

tree_file = sys.argv[1]
biom_file = sys.argv[2]
feature_dict_file = sys.argv[3]
sim_distance_file = sys.argv[4]

tree = dendropy.Tree.get_from_path(tree_file, schema='newick')
pdm = PhyloDM.load_from_dendropy(tree)
dm = pdm.dm(norm=False)
labels = pdm.taxa()

fid = biom.load_table(biom_file).ids(axis='observation')
prevalence = biom.load_table(biom_file).nonzero_counts(axis='observation')
feature_dict = {fid[i]: prevalence[i] for i in range(len(fid))}

with open(feature_dict_file, 'w') as f:
    for i in np.array(labels).astype(str):
        f.write(i + " " + str(feature_dict[i]) + '\n')

dt = np.dtype([('word1', 'int32'),
               ('word2', 'int32'),
               ('value', 'float64')])

data = coo_array(dm)
row = data.row
col = data.col
value = data.data
value = 1 - value / (np.max(value) + 1e-4)

temp = np.empty((len(row), ), dt)
for i in range(len(row)):
    temp[i] = (row[i] + 1, col[i] + 1, value[i])
temp.tofile(sim_distance_file)

