#!/bin/python

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier  
from sklearn.model_selection import cross_val_score 
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression

co_embedding = pd.read_csv("../../data/social_niche_embedding_100.txt",
                          header=None, sep=" ", low_memory=False, index_col=0)
co_embedding = co_embedding.drop("<unk>")

phy_embedding = pd.read_csv("../../data/Embedding_list/PCA_100.txt",
                          header=None, sep=" ", low_memory=False, index_col=0)
phy_embedding = phy_embedding.loc[co_embedding.index]


hgt_embed = pd.read_csv("Data/hgt.csv")

# 初始化随机森林分类器
rf_classifier = RandomForestClassifier(
    n_estimators=500,  # 树的数量
    max_depth=None,    # 树的最大深度（不限制）
    random_state=42,    # 随机种子，确保结果可复现
    n_jobs=-1
)

kf = KFold(n_splits=5, shuffle=True, random_state=42)
all_id = np.unique(list(hgt_embed.id_1.unique()) + list(hgt_embed.id_2.unique()))

test = []
group = []
labels = []
proba = []
for fold, (train_index, val_index) in enumerate(kf.split(all_id)):
    
    train_id = all_id[train_index]
    hgt_embed_train = hgt_embed.loc[[i in train_id for i in hgt_embed.id_1.values]]
    hgt_embed_train = hgt_embed_train.loc[[i in train_id for i in hgt_embed_train.id_2.values]]
    
    hgt_embed_test = hgt_embed.loc[[i  not in train_id for i in hgt_embed.id_1.values]]
    hgt_embed_test = hgt_embed_test.loc[[i not in train_id for i in hgt_embed_test.id_2.values]]

    hgt_embed_train_keep = hgt_embed_train.loc[hgt_embed_train.hgt > 0]
    hgt_embed_train_sample = hgt_embed_train.loc[hgt_embed_train.hgt == 0]

    hgt_embed_test_keep = hgt_embed_test.loc[hgt_embed_test.hgt > 0]
    hgt_embed_test_sample = hgt_embed_test.loc[hgt_embed_test.hgt == 0]

    data = pd.DataFrame({"id_1": hgt_embed_train_keep.id_1.tolist() + hgt_embed_train_sample.id_1.tolist(), 
                         "id_2": hgt_embed_train_keep.id_2.tolist() + hgt_embed_train_sample.id_2.tolist(),
                         "hgt": [1] * hgt_embed_train_keep.shape[0] + [0] * hgt_embed_train_sample.shape[0]})
    X_train = np.hstack((co_embedding.loc[data.id_1.values], co_embedding.loc[data.id_2.values]))
    y_train = data.hgt.values

    data = pd.DataFrame({"id_1": hgt_embed_test_keep.id_1.tolist() + hgt_embed_test_sample.id_1.tolist(), 
                         "id_2": hgt_embed_test_keep.id_2.tolist() + hgt_embed_test_sample.id_2.tolist(),
                         "hgt": [1] * hgt_embed_test_keep.shape[0] + [0] * hgt_embed_test_sample.shape[0]})
    x_test = np.hstack((co_embedding.loc[data.id_1.values], co_embedding.loc[data.id_2.values]))
    y_test = data.hgt.values

    # 训练模型
    rf_classifier.fit(X_train, y_train)
    # 对测试集进行预测
    y_pred_proba = rf_classifier.predict_proba(x_test)
    auc_score = roc_auc_score(y_test, y_pred_proba[:,1])
    
    test = test + [f"test_{fold}"] * len(y_test)
    group = group + ["SNE"] * len(y_test)
    labels = labels + list(y_test)
    proba = proba + list(y_pred_proba[:,1])
    
    print(f"test_{fold} AUC 值: {auc_score:.4f}")


    data = pd.DataFrame({"id_1": hgt_embed_train_keep.id_1.tolist() + hgt_embed_train_sample.id_1.tolist(), 
                     "id_2": hgt_embed_train_keep.id_2.tolist() + hgt_embed_train_sample.id_2.tolist(),
                     "hgt": [1] * hgt_embed_train_keep.shape[0] + [0] * hgt_embed_train_sample.shape[0]})
    X_train = np.hstack((phy_embedding.loc[data.id_1.values], phy_embedding.loc[data.id_2.values]))
    y_train = data.hgt.values

    data = pd.DataFrame({"id_1": hgt_embed_test_keep.id_1.tolist() + hgt_embed_test_sample.id_1.tolist(), 
                     "id_2": hgt_embed_test_keep.id_2.tolist() + hgt_embed_test_sample.id_2.tolist(),
                     "hgt": [1] * hgt_embed_test_keep.shape[0] + [0] * hgt_embed_test_sample.shape[0]})
    x_test = np.hstack((phy_embedding.loc[data.id_1.values], phy_embedding.loc[data.id_2.values]))
    y_test = data.hgt.values

    # 训练模型
    rf_classifier.fit(X_train, y_train) 
    # 对测试集进行预测
    y_pred_proba = rf_classifier.predict_proba(x_test)
    auc_score = roc_auc_score(y_test, y_pred_proba[:,1])
    print(f"AUC_{fold} 值: {auc_score:.4f}")

    test = test + [f"test_{fold}"] * len(y_test)
    group = group + ["PhyloE"] * len(y_test)
    labels = labels + list(y_test)
    proba = proba + list(y_pred_proba[:,1])


hgt_predict_res = pd.DataFrame({"test":test, "group":group, "labels":labels, "proba":proba})
hgt_predict_res.to_csv("Data/hgt_predict_res_all.csv", index=None)

