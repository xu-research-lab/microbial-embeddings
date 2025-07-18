import sys
import shap
import biom
import torch
import random
from torch import nn
import torch.nn.functional as F
import collections
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

from torch.autograd import Variable
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import confusion_matrix, f1_score, matthews_corrcoef, mean_absolute_error, r2_score
from sklearn import metrics

from membed.Attention_embedding import Fid, load_data_imdb, DataLoader, TransformerEncoder


def load_array(data_arrays, batch_size):
    """Construct a PyTorch data iterator.
    Defined in :numref:`sec_utils`"""
    dataset = torch.utils.data.TensorDataset(*data_arrays)
    return torch.utils.data.DataLoader(dataset, batch_size, shuffle=False)
    

def load_data_2(otu):

    features = np.zeros((otu.shape[0], num_steps), dtype=int)
    cls_column = np.zeros((features.shape[0], 1)).astype(int)
    abundance = np.zeros((otu.shape[0], num_steps))
    cls_abundance = np.ones((features.shape[0], 1)).astype(int)
    
    for i in range(0, otu.shape[0]):
        nonzero_count = np.count_nonzero(otu[i,])
        if nonzero_count >= num_steps:
            sorted_indices = np.argsort(otu[i,])[::-1]
            features[i, ] = np.array([fid_dict[line]
                                    for line in
                                    fid[sorted_indices[:num_steps]]])
            abundance[i, ] = otu[i, sorted_indices[:num_steps]]
        else:
            indices = np.nonzero(otu[i,])
            features[i, :nonzero_count] = np.array([fid_dict[line]
                                                    for line in
                                                    fid[indices]])
            features[i, nonzero_count:] = fid_dict['<pad>']
            abundance[i, :nonzero_count] = otu[i, indices]


    features = np.concatenate((cls_column, features), axis=1)
    abundance = np.concatenate((cls_abundance, abundance), axis=1)

    input_mask = np.ones((features.shape[0], features.shape[1]))
    temp = features == fid_dict['<pad>']
    input_mask[temp] = 0 # mask padding tokens
    input_mask[features == 0] = 0 # mask cls tokens
    
    return load_array((torch.tensor(features), torch.tensor(abundance), torch.tensor(input_mask)), 900)


num_steps = 600
batch_size = 6000
d_model = 100
n_layers = 1
n_heads = 1
p_drop = 0.4
d_ff = 8
group = "group"

file_path = "../../Data/disease_data/IBD"
IBD_study = ["PRJNA324147", "PRJNA368966", "PRJNA422193", "PRJNA431126", "PRJNA450340", "qiita_1629", "qiita_2538", "RISK_PRISM_f"]
metadata = "../../Data/disease_data/IBD/metadata.tsv"

i = int(sys.argv[1])
cuda_id = sys.argv[2]
devices = torch.device(f'cuda:{cuda_id}')

train_table = biom.load_table(f"{file_path}/train_{i+1}.biom")
train_table = train_table.rankdata(axis='sample', inplace=False)
table_1 = train_table.matrix_data.multiply(1 / train_table.max(axis='sample'))
table_1 = table_1.toarray().T
table_1 = pd.DataFrame(data=table_1, index=train_table.ids(axis='sample'),
                       columns=train_table.ids(axis='observation'))

X_valid_table = biom.load_table(f"{file_path}/test_{i+1}.biom")
X_valid_table = X_valid_table.rankdata(axis='sample', inplace=False)
X_valid = X_valid_table.matrix_data.multiply(1 / X_valid_table.max(axis='sample'))
X_valid = X_valid.toarray().T
X_valid = pd.DataFrame(data=X_valid, index=X_valid_table.ids(axis='sample'),
                       columns=X_valid_table.ids(axis='observation'))

table = pd.concat([table_1, X_valid])
# table = table_1
med = pd.DataFrame(data=table.median().values.reshape((1, table.shape[1])),
                   columns=table.columns)

# 加载模型
train_data = load_data_imdb(f"{file_path}/train_{i+1}.biom", metadata, group, num_steps)
train_iter = DataLoader(train_data, batch_size=batch_size, shuffle=False)

fid_dict = train_data()
net = TransformerEncoder(otu_size=len(fid_dict),
                         seq_len=num_steps+1,
                         d_model=d_model,
                         n_layers=n_layers,
                         n_heads=n_heads,
                         p_drop=p_drop,
                         d_ff=d_ff,
                         pad_id=fid_dict['<pad>'])

net.load_state_dict( 
    torch.load(
        f"{file_path}/attention_{i+1}.pt",
        map_location=torch.device('cpu')))

fid = table.columns

def f(x):
    """
    x:numpy array
    """
    data_iter = load_data_2(np.array(x))
    net_copy = net
    net_copy = net_copy.to(devices)
    with torch.no_grad():
        net_copy.eval()
        
        for i, (features, abundance, mask) in enumerate(data_iter):
            features = features.to(devices)
            abundance = abundance.to(devices)
            mask = mask.to(devices)
            pred, _ = net_copy(features.int(), abundance, mask)
            if i == 0:
                prob = pred
            else:
                prob = torch.cat((prob, pred))
    prob = torch.squeeze(prob, dim=1)
    del features, abundance, pred, net_copy
    prob = prob.cpu().detach().numpy()

    return prob


explainer = shap.Explainer(f, med)
shap_values = explainer(table, max_evals=9000)
shap_table = pd.DataFrame(data=shap_values.values, index=table.index, columns=table.columns)

shap_table.to_csv(f"../../Data/biomark/shap_table_{IBD_study[i]}.csv", index=True)

