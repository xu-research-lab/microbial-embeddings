#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
from phylodm import PhyloDM
from sklearn.decomposition import PCA
# import umap.umap_ as umap
import os


def main(tree_path, n_dims, embedding_file):
    # read tree file, format:newick
    pdm = PhyloDM.load_from_newick_path(tree_path)
    # obtain pairwise distance matrix and labels
    mat = pdm.dm(norm=False)
    labels = pdm.taxa()
    # reduce the memory possession
    mat = mat.astype("float16")
    # reduce the memory possession
    del pdm
    # pca reduction
    pca = PCA(n_components=n_dims) # 指定降维后的维度
    reduced_pca = pca.fit_transform(mat)
    reduced_pca = pd.DataFrame(reduced_pca, index=labels)
    new_row = pd.Series(np.zeros(reduced_pca.shape[1]), name='<unk>')
    reduced_pca = reduced_pca.append(new_row)
    reduced_pca.round(4).to_csv(embedding_file, sep=" ", header=False)
    # np.savetxt(os.path.join('/home/wangbin/phy2vec/data/tradition', tree_name + '_' + 'embeddings' + '_' + 'pca' + '.txt'), 
    #            reduced_pca, delimiter="\t", fmt="%.4f")
    # del reduced_pca, pca
    # # umap reduction
    # umaper = umap.UMAP(n_components=n_dims)  # 指定降维后的维度
    # reduced_umap = umaper.fit_transform(mat) 
    # np.savetxt(os.path.join('/home/wangbin/phy2vec/data/tradition',tree_name + '_' + 'embeddings' + '_' + 'umap' + '.txt'), reduced_umap, delimiter="\t", fmt="%.4f")
    # del reduced_umap, umaper

    # np.savetxt(os.path.join('/home/wangbin/phy2vec/data/tradition', tree_name + '_' + 'labels' + '_' + 'pca_uamp' + '.txt'), labels, delimiter="\t", fmt="%s")

if __name__ == '__main__':
    tree_file = sys.argvs[1]
    n_dims = sys.argvs[2]
    embedding_file = sys.argvs[3]
    main(tree_file, n_dims, embedding_file)
