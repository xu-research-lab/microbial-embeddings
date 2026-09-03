#!/bin/python

import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier  
from sklearn.model_selection import cross_val_score 
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression


def read_cooccur(cooccur_file, feature_dict):
    """Load GloVe co-occurrence records as a {(fid_a, fid_b): value} dict.

    cooccur_file: headerless binary array of dtype
        [('otu1','i4'),('otu2','i4'),('value','f8')], 1-based feature indices,
        symmetrically expanded. See cooccur_embedding.cooccur_workflow.
    feature_dict: space-separated "<feature_id> <count>" lines in GloVe vocab
        order; line N is feature index N (1-based). See
        cooccur_embedding.get_feature_dict.
    """
    fids = pd.read_csv(feature_dict, sep=" ", header=None)[0].tolist()
    idx2fid = {i + 1: f for i, f in enumerate(fids)}

    dt = np.dtype([('otu1', 'i4'), ('otu2', 'i4'), ('value', 'f8')])
    rec = np.fromfile(cooccur_file, dtype=dt)
    return {(idx2fid[a], idx2fid[b]): v
            for a, b, v in zip(rec['otu1'], rec['otu2'], rec['value'])}


def cooccur_feature(cooccur, id_1, id_2):
    """Column vector of the pair co-occurrence value (0 when the pair is absent)."""
    return np.array([[cooccur.get((a, b), cooccur.get((b, a), 0.0))]
                     for a, b in zip(id_1, id_2)])


# TODO: real paths -- .cooccur binary from cooccur_embedding.cooccur_workflow
# and the vocab file from cooccur_embedding.get_feature_dict.
COOCCUR_FILE = "../../script/cooccur_otuembedding/table.co"
FEATURE_DICT = "../../script/cooccur_otuembedding/feature-dict.csv"

co_embedding = pd.read_csv("../../data/social_niche_embedding_100.txt",
                          header=None, sep=" ", low_memory=False, index_col=0)
co_embedding = co_embedding.drop("<unk>")

phy_embedding = pd.read_csv("../../data/phylo_embed_PCA_100.txt",
                          header=None, sep=" ", low_memory=False, index_col=0)
phy_embedding = phy_embedding.loc[co_embedding.index]


hgt_embed = pd.read_csv("data/hgt.csv")

cooccur = read_cooccur(COOCCUR_FILE, FEATURE_DICT)

rf_classifier = RandomForestClassifier(
    n_estimators=500,  
    max_depth=None,    
    random_state=42,  
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

    rf_classifier.fit(X_train, y_train)
    y_pred_proba = rf_classifier.predict_proba(x_test)
    auc_score = roc_auc_score(y_test, y_pred_proba[:,1])
    
    test = test + [f"test_{fold}"] * len(y_test)
    group = group + ["SNE"] * len(y_test)
    labels = labels + list(y_test)
    proba = proba + list(y_pred_proba[:,1])
    
    print(f"test_{fold} AUC: {auc_score:.4f}")


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

    rf_classifier.fit(X_train, y_train) 
    y_pred_proba = rf_classifier.predict_proba(x_test)
    auc_score = roc_auc_score(y_test, y_pred_proba[:,1])
    print(f"AUC_{fold}: {auc_score:.4f}")

    test = test + [f"test_{fold}"] * len(y_test)
    group = group + ["PhyloE"] * len(y_test)
    labels = labels + list(y_test)
    proba = proba + list(y_pred_proba[:,1])


    data = pd.DataFrame({"id_1": hgt_embed_train_keep.id_1.tolist() + hgt_embed_train_sample.id_1.tolist(),
                     "id_2": hgt_embed_train_keep.id_2.tolist() + hgt_embed_train_sample.id_2.tolist(),
                     "hgt": [1] * hgt_embed_train_keep.shape[0] + [0] * hgt_embed_train_sample.shape[0]})
    X_train = cooccur_feature(cooccur, data.id_1.values, data.id_2.values)
    y_train = data.hgt.values

    data = pd.DataFrame({"id_1": hgt_embed_test_keep.id_1.tolist() + hgt_embed_test_sample.id_1.tolist(),
                     "id_2": hgt_embed_test_keep.id_2.tolist() + hgt_embed_test_sample.id_2.tolist(),
                     "hgt": [1] * hgt_embed_test_keep.shape[0] + [0] * hgt_embed_test_sample.shape[0]})
    x_test = cooccur_feature(cooccur, data.id_1.values, data.id_2.values)
    y_test = data.hgt.values

    rf_classifier.fit(X_train, y_train)
    y_pred_proba = rf_classifier.predict_proba(x_test)
    auc_score = roc_auc_score(y_test, y_pred_proba[:,1])
    print(f"cooccur_{fold} AUC: {auc_score:.4f}")

    test = test + [f"test_{fold}"] * len(y_test)
    group = group + ["Cooccur"] * len(y_test)
    labels = labels + list(y_test)
    proba = proba + list(y_pred_proba[:,1])


    data_train = pd.DataFrame({"id_1": hgt_embed_train_keep.id_1.tolist() + hgt_embed_train_sample.id_1.tolist(),
                     "id_2": hgt_embed_train_keep.id_2.tolist() + hgt_embed_train_sample.id_2.tolist(),
                     "hgt": [1] * hgt_embed_train_keep.shape[0] + [0] * hgt_embed_train_sample.shape[0]})
    train_emb_id1 = pd.concat([co_embedding.loc[data_train.id_1.values], phy_embedding.loc[data_train.id_1.values]], axis=1)
    train_emb_id2 = pd.concat([co_embedding.loc[data_train.id_2.values], phy_embedding.loc[data_train.id_2.values]], axis=1)
    X_train = np.hstack([train_emb_id1.values, train_emb_id2.values])
    y_train = data_train.hgt.values

    data_test = pd.DataFrame({"id_1": hgt_embed_test_keep.id_1.tolist() + hgt_embed_test_sample.id_1.tolist(),
                     "id_2": hgt_embed_test_keep.id_2.tolist() + hgt_embed_test_sample.id_2.tolist(),
                     "hgt": [1] * hgt_embed_test_keep.shape[0] + [0] * hgt_embed_test_sample.shape[0]})
    test_emb_id1 = pd.concat([co_embedding.loc[data_test.id_1.values], phy_embedding.loc[data_test.id_1.values]], axis=1)
    test_emb_id2 = pd.concat([co_embedding.loc[data_test.id_2.values], phy_embedding.loc[data_test.id_2.values]], axis=1)
    x_test = np.hstack([test_emb_id1.values, test_emb_id2.values])
    y_test = data_test.hgt.values

    rf_classifier.fit(X_train, y_train)
    y_pred_proba = rf_classifier.predict_proba(x_test)
    auc_score = roc_auc_score(y_test, y_pred_proba[:,1])
    print(f"SNE+PhyloE_{fold} AUC: {auc_score:.4f}")

    test = test + [f"test_{fold}"] * len(y_test)
    group = group + ["SNE+PhyloE"] * len(y_test)
    labels = labels + list(y_test)
    proba = proba + list(y_pred_proba[:,1])


hgt_predict_res = pd.DataFrame({"test":test, "group":group, "labels":labels, "proba":proba})
hgt_predict_res.to_csv("data/hgt_predict_res_all.csv", index=None)
