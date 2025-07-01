#!/bin/python

import sys
import biom
import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.metrics import matthews_corrcoef
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score, f1_score, fbeta_score, accuracy_score
from sklearn.model_selection import StratifiedShuffleSplit

from sklearn.metrics import confusion_matrix
from inspect import signature

train_file = sys.argv[1]
test_file = sys.argv[2]
metadata = sys.argv[3]
embedding_file = sys.argv[4]
group = sys.argv[5]
plot_file = sys.argv[6]
ROC_file = sys.argv[7]
PREDICTED_PROBS_FILE = sys.argv[8]

def computeMLstats(m,
                   data,
                   y,
                   sample_ids=None,             # 新增：样本ID参数
                   probs_output_file=None,    # 新增：概率输出文件名参数
                   plot=False,
                   plot_file=False,
                   plot_pr=False,
                   graph_title=None,
                   flipped=False):
        
    probs = m.predict_proba(data)
    
    if sample_ids is not None and probs_output_file is not None:
        if len(sample_ids) != len(y):
            raise ValueError("样本ID的数量必须与真实标签y的数量一致。")
        if len(sample_ids) != probs.shape[0]:
            raise ValueError("样本ID的数量必须与输入数据data的样本数量一致。")

    # 创建DataFrame保存结果
    # probs[:, 1] 是模型预测样本属于类别1的概率
    # probs[:, 0] 是模型预测样本属于类别0的概率 (如果需要也可以保存)
    predictions_df = pd.DataFrame({
        'SampleID': sample_ids,
        'TrueLabel': y,
        'PredictedProbability': probs[:, 1]
        # 如果需要类别0的概率: 'PredictedProbability_Class0': probs[:, 0]
    })
    predictions_df.to_csv(probs_output_file, index=False)
    print(f"样本预测概率已保存到: {probs_output_file}")

    # Flip for opposite class imbalance
    if flipped:
        y = [1 - i for i in y]
        probs = 1 - probs
    # Compute ROC curve and area the curve
    fpr, tpr, thresholds = roc_curve(y, probs[:, 1], pos_label=1)
    roc_auc = auc(fpr, tpr)  # 等价于 roc_auc_score(truth,prob)
    # Compute precision-recall
    precision, recall, threshold = precision_recall_curve(y,
                                                          probs[:, 1],
                                                          pos_label=1)

    average_precision = average_precision_score(y, probs[:, 1])
    numerator = 2 * recall * precision
    denom = recall + precision
    f1_scores = np.divide(numerator,
                          denom,
                          out=np.zeros_like(denom),
                          where=(denom != 0))

    max_f1 = np.max(f1_scores)
    max_f1_thresh = threshold[np.argmax(f1_scores)]
    aupr = metrics.auc(recall, precision)
    # max_f1_thresh = 0.5
    y_hat = np.array([1 if i > max_f1_thresh else 0 for i in probs[:, 1]])
    mcc = matthews_corrcoef(y, y_hat)
    cm = confusion_matrix(y, y_hat)
    f1 = f1_score(y, y_hat, average='macro')
    acc = accuracy_score(y, y_hat)
    f2 = fbeta_score(y, np.argmax(probs, axis=1), beta=2)
    if plot:
        plt.subplot(1, 2, 1)
        plt.plot(fpr, tpr, lw=2, alpha=0.3, label='AUC ROC = %0.2f' % roc_auc)
        # 'AUC PR = %0.2f' % pr_avg_pr
        plt.legend(loc="lower right")
        x = np.linspace(0, 1, 10)
        plt.plot(x, x)
        plt.title(graph_title)
        plt.xlabel('False positive rate')
        plt.ylabel('True positive rate')
        plt.savefig(plot_file)

    if plot_pr:
        plt.subplot(1, 2, 2)
        # In matplotlib < 1.5, plt.fill_between does not have a 'step' argument
        step_kwargs = ({
            'step': 'post'
        } if 'step' in signature(plt.fill_between).parameters else {})
        plt.step(recall,
                 precision,
                 color='b',
                 alpha=0.2,
                 where='post',
                 label='AUC PR = %0.2f' % average_precision)
        plt.fill_between(recall,
                         precision,
                         alpha=0.2,
                         color='b',
                         **step_kwargs)
        plt.legend(loc="lower right")
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
    return roc_auc, fpr, tpr, average_precision, f1, f2, y_hat, mcc, aupr, acc


def getFeatureImportance(m, data, y):
    feat_imp = m.feature_importances_
    feat_imp_labeled = zip(data.columns.values, feat_imp)
    feat_imp_sort = sorted(feat_imp_labeled, key=lambda t: t[1], reverse=True)
    return feat_imp_sort


def predictIBD(X_train,
               y_train,
               X_test,
               y_test,
               graph_title="",
               max_depth=12,
               n_estimators=140,
               plot=False,
               plot_pr=False,
               weight=20,
               feat_imp=False,
               flipped=False):
    global feat_imp_sort
    weights = {0: 1, 1: weight}  ## ？
    m = RandomForestClassifier(max_depth=max_depth,
                               random_state=0,
                               n_estimators=n_estimators,
                               class_weight=weights,
                               n_jobs=24)
    m.fit(X_train, y_train)
    # fpr: false positive rates
    # tpr: true positive rates
    roc_auc, fpr, tpr, average_precision, f1, f2 = computeMLstats(
        m,
        data=X_test,
        y=y_test,
        plot=plot,
        plot_pr=plot_pr,
        graph_title=graph_title,
        flipped=flipped)
    if feat_imp:
        feat_imp_sort = getFeatureImportance(m, data=X_train, y=y_train)

    return m, roc_auc, fpr, tpr, average_precision, f1, f2, feat_imp_sort


def crossValPrediction(otu_use,
                       y,
                       max_depth=10,
                       n_estimators=65,
                       weight=5,
                       plot=True,
                       plot_pr=True,
                       folds=5):
    kf = StratifiedShuffleSplit(n_splits=folds)
    kf.get_n_splits(otu_use, y)

    AU_ROC, AP = [], []
    f1_cross_val, f2_cross_val = [], []
    feat_imp_crossVal = []
    models = []
    i = 0
    for train_index, val_index in kf.split(otu_use, y):
        otu_train = otu_use.iloc[train_index, :]
        otu_val = otu_use.iloc[val_index, :]
        y_train = np.array(y)[train_index]
        y_val = np.array(y)[val_index]

        if plot and plot_pr:
            plt.subplot(1, 2, 1)
        model, AUROC, fpr, tpr, average_prec, f1, f2, feat_imp = predictIBD(
            otu_train,
            y_train,
            otu_val,
            y_val,
            max_depth=max_depth,
            n_estimators=n_estimators,
            weight=weight,
            plot=plot,
            plot_pr=plot_pr,
            feat_imp=True)
        models.append(model)
        AU_ROC.append(AUROC)
        AP.append(average_prec)
        f1_cross_val.append(f1)
        f2_cross_val.append(f2)
        feat_imp_crossVal.append(feat_imp)
        i = i + 1
    return models, AU_ROC, AP, f1_cross_val, f2_cross_val, feat_imp_crossVal


def trainHyperParameters(X_train, y_train):
    depths = [2, 3, 5, 7, 10]
    n_estimators = [50, 65, 80, 95, 110, 125, 140, 155]
    weights = [1, 3, 5, 10, 15]

    df = np.zeros((len(depths) * len(n_estimators) * len(weights), 6))
    i = 0
    for depth in depths:
        for trees in n_estimators:
            for weight in weights:
                # models, AU_ROC, AP, f1_cross_val, f2_cross_val, feat_imp_crossVal
                _, auc, average_prec, f1_value, _, _ = crossValPrediction(
                    otu_use=X_train,
                    y=y_train,
                    max_depth=depth,
                    n_estimators=trees,
                    weight=weight,
                    plot=False,
                    plot_pr=False,
                    folds=5)
                df[i, :] = [
                    np.mean(auc),
                    np.mean(average_prec),
                    np.mean(f1_value), depth, trees, weight
                ]
                i += 1
    return df


train = biom.load_table(train_file)
test = biom.load_table(test_file)

train_col = train.ids(axis='observation')
test_col = test.ids(axis='observation')
train_indx = train.ids(axis='sample')
test_indx = test.ids(axis='sample')

train = train.rankdata(axis='sample', inplace=False)
train = train.matrix_data.multiply(1 / train.max(axis='sample'))
test = test.rankdata(axis='sample', inplace=False)
test = test.matrix_data.multiply(1 / test.max(axis='sample'))

# train = train.norm(axis='sample', inplace=False).matrix_data
# test = test.norm(axis='sample', inplace=False).matrix_data

otu_train = pd.DataFrame(columns=train_col,
                         index=train_indx,
                         data=train.toarray().T)
otu_test = pd.DataFrame(columns=test_col,
                        index=test_indx,
                        data=test.toarray().T)
map_keep = pd.read_csv(metadata, sep="\t", index_col=0, low_memory=False)

otu_train = otu_train.loc[
    np.intersect1d(otu_train.index.values, map_keep.index.values), ]

otu_test = otu_test.loc[
    np.intersect1d(otu_test.index.values, map_keep.index.values),
    otu_train.columns]
if embedding_file != 'None':
    embedding_vector = pd.read_csv(embedding_file,
                                   header=None,
                                   sep=" ",
                                   index_col=0,
                                   low_memory=False)
    embedding_vector.index = [str(i) for i in embedding_vector.index]
    fid = np.intersect1d(otu_train.columns.values,
                         embedding_vector.index.values)
    otu_train = otu_train.loc[:, fid].dot(embedding_vector.loc[fid, ])
    otu_test = otu_test.loc[:, fid].dot(embedding_vector.loc[fid, ])

X_train = np.array(otu_train)
x_test = np.array(otu_test)
y_train = np.array(map_keep.loc[np.intersect1d(otu_train.index.values,
                                               map_keep.index.values)][group])
y_test = np.array(map_keep.loc[np.intersect1d(otu_test.index.values,
                                              map_keep.index.values)][group])
test_sample_ids = np.intersect1d(otu_test.index.values, map_keep.index.values)

# df = trainHyperParameters(pd.DataFrame(X_train), y_train)
# df = pd.DataFrame(df,
#                   columns=['AUROC', 'AP', 'F1', 'depth', 'n_trees', 'weight'])
# df = df.sort_values(by=['AUROC', 'AP', 'F1'], ascending=False)

# area under the receiver operating curve (AUROC), area under the precision-recall curve (AUPR)
# 选取最优参数的最优模型
# models, AU_ROCs, APs, f1_values, f2_values, feat_imp_lists = crossValPrediction(
#     otu_use=pd.DataFrame(X_train),
#     y=y_train,
#     max_depth=int(df['depth'].values[0]),
#     n_estimators=int(df['n_trees'].values[0]),
#     weight=int(df['weight'].values[0]),
#     plot=False,
#     plot_pr=False,
#     folds=5)
# model = models[AU_ROCs.index(max(AU_ROCs))]
# roc_auc, fpr, tpr, average_precision, f1, f2
model = RandomForestClassifier(n_estimators=500,
                               n_jobs=24)
model.fit(pd.DataFrame(X_train), y_train)
AUC_ROC, fpr, tpr, average_precision, f1, f2, y_test_pred, mcc, aupr, acc = computeMLstats(
    model,
    data=pd.DataFrame(x_test),
    y=y_test,
    sample_ids=test_sample_ids,    # 传入样本ID
    probs_output_file=PREDICTED_PROBS_FILE, # 传入概率输出文件名
    plot=True,
    plot_file=plot_file,
    plot_pr=False,
    flipped=False)
print("The AUC ROC: %0.4f, The Average Precision: %0.4f, The F1 value: %0.4f, The mcc value: %0.4f, The aupr value: %0.4f, The acc value: %0.4f" %
      (AUC_ROC, average_precision, f1, mcc, aupr, acc))

roc_data = pd.DataFrame({
    "FPR": fpr,
    "TPR": tpr
})
roc_data.to_csv(ROC_file, index=False)


# 此时的混淆矩阵是二分类的
Conf_mat = confusion_matrix(y_test, list(y_test_pred), labels=[0.0, 1.0])
print("The conflusion Matrix is:\n", Conf_mat)
print(metrics.classification_report(y_test_pred, y_test))
