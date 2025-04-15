import math
import time
import torch
import biom
import collections
import numpy as np
import pandas as pd
from torch import nn
import torch.nn.functional as F
from torch.autograd import Variable
from matplotlib import pyplot as plt
from sklearn.metrics import confusion_matrix, f1_score, matthews_corrcoef, mean_absolute_error, r2_score
from sklearn import metrics
from sklearn.utils.class_weight import compute_class_weight

import gc
import os

import yaml


def read_imdb(biom_table, embedding_file, metadata, group):
    """读取OTU table和样本标签"""
    labels = []
    table = biom.load_table(biom_table)
    fid = table.ids(axis="observation")
    if embedding_file is not None:
        embedding = pd.read_csv(embedding_file,
                                index_col=0,
                                low_memory=False,
                                sep=" ",
                                header=None)
        embedding.index = [str(i) for i in embedding.index]
        id = np.intersect1d(fid, embedding.index.values)

        no_embedding_id = fid[[False if i in id else True for i in fid]]
        for i in no_embedding_id:
            embedding.loc[i] = np.random.uniform(-0.5, 0.5, embedding.shape[1])
        embedding = embedding.loc[fid, ]
    # table = table.filter(id, axis='observation', inplace=False)
    table = table.rankdata(axis='sample', inplace=False)
    abs_ = table.matrix_data.multiply(1 / table.max(axis='sample'))
    if embedding_file is not None:
        otu = pd.DataFrame(
            columns=table.ids(axis='observation'),
            index=table.ids(axis='sample'),
            data=abs_.toarray().T)
        otu = np.array(otu.dot(embedding))
    else:
        otu = abs_.toarray().T
    mapping_file = pd.read_csv(metadata,
                               sep="\t",
                               index_col=0,
                               low_memory=False)                     
    labels = list(mapping_file.loc[table.ids(axis='sample')][group].values)
    return otu, labels


def load_array(data_arrays, batch_size, is_train=True):
    """Construct a PyTorch data iterator.
    Defined in :numref:`sec_utils`"""
    dataset = torch.utils.data.TensorDataset(*data_arrays)
    return torch.utils.data.DataLoader(dataset, batch_size, shuffle=is_train)


def load_data_imdb(batch_size, train_biom, test_biom, metadata, embedding_file,
                   group):
    """Return data iterators and the vocabulary of the IMDb review dataset.
    Defined in :numref:`sec_sentiment`

    Returns
    -------
    _type_
        _description_ 
    """
    otu, labels = read_imdb(train_biom, embedding_file, metadata, group)
    train_data = (otu, labels)
    test_data = read_imdb(test_biom, embedding_file, metadata, group)
    train_iter = load_array(
        (torch.tensor(train_data[0]), torch.tensor(train_data[1])), batch_size)
    test_iter = load_array(
        (torch.tensor(test_data[0]), torch.tensor(test_data[1])),
        batch_size,
        is_train=False)
    return train_iter, test_iter, labels


def cpu():
    """Defined in :numref:`sec_use_gpu`"""
    return torch.device('cpu')


def gpu(i=0):
    """Defined in :numref:`sec_use_gpu`"""
    return torch.device(f'cuda:{i}')


def num_gpus():
    """Defined in :numref:`sec_use_gpu`"""
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = "1"
    return torch.cuda.device_count()


def try_all_gpus(numb):
    """Return all available GPUs, or [cpu(),] if no GPU exists.
    Defined in :numref:`sec_use_gpu`"""
    return [gpu(i) for i in [numb]]


def try_gpu(i=0):
    """Return gpu(i) if exists, otherwise return cpu().
    Defined in :numref:`sec_use_gpu`"""
    if num_gpus() >= i + 1:
        return gpu(i)
    return cpu()


def evaluate_auc_gpu(Y, prob, device=None):
    """_summary_

    Parameters
    ----------
    net : _type_
        _description_
    data_iter : _type_
        _description_
    device : _type_, optional
        _description_, by default None
    """
    Y = Y.numpy().astype('int')
    prob = prob.cpu().detach().numpy()
    fpr, tpr, thresholds = metrics.roc_curve(Y, prob, pos_label=1)
    auc = metrics.auc(fpr, tpr)
    precision, recall, threshold = metrics.precision_recall_curve(Y,
                                                                  prob,
                                                                  pos_label=1)
    numerator = 2 * recall * precision
    denom = recall + precision

    f1_scores = np.divide(numerator,
                          denom,
                          out=np.zeros_like(denom),
                          where=(denom != 0))
    max_f1 = np.max(f1_scores)
    max_f1_thresh = threshold[np.argmax(f1_scores)]
    # max_f1_thresh = 0.5
    aupr = metrics.auc(recall, precision)
    y_hat = np.array([1 if i > max_f1_thresh else 0 for i in prob])
    cm = confusion_matrix(Y, y_hat)
    f1_scores = f1_score(Y, y_hat, average='macro')
    mcc = matthews_corrcoef(Y, y_hat)

    return auc, max_f1_thresh, aupr, cm, f1_scores, mcc


def evaluate_r2_gpu(Y, pred,device=None):
    Y = Y.numpy() #.astype('int')
    pred = pred.cpu().detach().numpy()
    mae = mean_absolute_error(Y, pred)
    r2 = r2_score(Y, pred)

    return mae, r2


class Timer:
    """Record multiple running times."""

    def __init__(self):
        """Defined in :numref:`sec_minibatch_sgd`"""
        self.times = []
        self.start()

    def start(self):
        """Start the timer."""
        self.tik = time.time()

    def stop(self):
        """Stop the timer and record the time in a list."""
        self.times.append(time.time() - self.tik)
        return self.times[-1]

    def avg(self):
        """Return the average time."""
        return sum(self.times) / len(self.times)

    def sum(self):
        """Return the sum of time."""
        return sum(self.times)

    def cumsum(self):
        """Return the accumulated time."""
        return np.array(self.times).cumsum().tolist()


def train_batch_ch13(net, X, y, loss, trainer, devices):
    """Train for a minibatch with multiple GPUs (defined in Chapter 13).
    Defined in :numref:`sec_image_augmentation`"""
    if isinstance(X, list):
        # Required for BERT fine-tuning (to be covered later)
        X = [x.to(devices[0]) for x in X]
    else:
        X = X.to(devices[0])
    y = y.to(devices[0])
    net.train()
    trainer.zero_grad()
    pred = net(X.float())
    pred = torch.squeeze(pred, dim=1)
    l = loss(pred, y.float())
    l.backward()
    trainer.step()
    train_loss_sum = l.sum()
    return train_loss_sum, pred


class Net(nn.Module):

    def __init__(self, feature_size, num_hiddens, p_drop, **kwargs):
        super(Net, self).__init__()
        self.fc1 = nn.Linear(feature_size, num_hiddens)
        self.fc2 = nn.Linear(num_hiddens, num_hiddens)
        self.fc3 = nn.Linear(num_hiddens, 1)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()
        self.dropout = nn.Dropout(p=p_drop)

    def forward(self, x, classification=True):
        x = self.relu(self.fc1(x))
        x = self.dropout(x)
        # x = self.relu(self.fc2(x))
        # x = self.dropout(x)
        if classification is True:
            x = self.sigmoid(self.fc3(x))
        else:
            x = self.relu(self.fc3(x))

        return x
    

def train_log(epoch, **kwargs):
    wandb.log(kwargs)


def set_axes(axes, xlabel, ylabel, xlim, ylim, xscale, yscale, legend):
    """Set the axes for matplotlib.

    Defined in :numref:`sec_calculus`"""
    axes.set_xlabel(xlabel), axes.set_ylabel(ylabel)
    axes.set_xscale(xscale), axes.set_yscale(yscale)
    axes.set_xlim(xlim), axes.set_ylim(ylim)
    if legend:
        axes.legend(legend)
    axes.grid()


class Animator:
    """For plotting data in animation."""

    def __init__(self,
                 xlabel=None,
                 ylabel=None,
                 legend=None,
                 xlim=None,
                 ylim=None,
                 xscale='linear',
                 yscale='linear',
                 fmts=('-', 'm--', 'g-.', 'r:'),
                 nrows=1,
                 ncols=1,
                 figsize=(3.5, 2.5)):
        """Defined in :numref:`sec_utils`"""
        # Incrementally plot multiple lines
        if legend is None:
            legend = []
        # use_svg_display()
        self.fig, self.axes = plt.subplots(nrows, ncols, figsize=figsize)
        if nrows * ncols == 1:
            self.axes = [
                self.axes,
            ]
        # Use a lambda function to capture arguments
        self.config_axes = lambda: set_axes(self.axes[0], xlabel, ylabel, xlim,
                                            ylim, xscale, yscale, legend)
        self.X, self.Y, self.fmts = None, None, fmts

    def add(self, x, y, plotfile):
        # Add multiple data points into the figure
        if not hasattr(y, "__len__"):
            y = [y]
        n = len(y)
        if not hasattr(x, "__len__"):
            x = [x] * n
        if not self.X:
            self.X = [[] for _ in range(n)]
        if not self.Y:
            self.Y = [[] for _ in range(n)]
        for i, (a, b) in enumerate(zip(x, y)):
            if a is not None and b is not None:
                self.X[i].append(a)
                self.Y[i].append(b)
        self.axes[0].cla()
        for x, y, fmt in zip(self.X, self.Y, self.fmts):
            self.axes[0].plot(x, y, fmt)
        self.config_axes()
        self.fig.savefig(plotfile)


def train_cls(net, train_iter, test_iter, loss, trainer, num_epochs, devices, 
              plotfile_loss, plotfile_auc, mlp_model):
    animator_1 = Animator(xlabel='epoch', xlim=[1, num_epochs],
                          legend=['train loss', 'test loss'])
    animator_2 = Animator(xlabel='epoch', xlim=[1, num_epochs], ylim=[0, 1],
                          legend=['train auc', 'test auc'])
    test_list = []
    train_list = []
    test_loss_list = []
    for epoch in range(num_epochs):
        train_loss = 0
        for i, (features, labels) in enumerate(train_iter):
            l, pred = train_batch_ch13(net, features, labels, loss, trainer, devices)
            train_loss += l
            if i == 0:
                prob1 = pred
                Label = labels
            else:
                prob1 = torch.cat((prob1, pred), dim=0)
                Label = torch.cat((Label, labels), dim=0)
        train_loss = train_loss / (i+1)
        train_auc, train_f1_thresh, train_aupr, train_cm, train_f1, train_mcc = evaluate_auc_gpu(
            Label, prob1)
        animator_1.add(epoch + 1, (train_loss.cpu().detach().numpy(), None), plotfile_loss)
        animator_2.add(epoch + 1, (train_auc, None), plotfile_auc)
        with torch.no_grad():
            net.eval()
            test_loss = 0
            for i, (features, labels) in enumerate(test_iter):
                X = features.to(devices[0])
                pred = net(X.float())
                y = labels.to(devices[0])
                test_loss += loss(torch.squeeze(pred, dim=1), y.float())
                if i == 0:
                    prob2 = pred
                    Label = labels
                else:
                    prob2 = torch.cat((prob2, pred), dim=0)
                    Label = torch.cat((Label, labels), dim=0)
            test_loss = test_loss / (i +1)
        test_auc, test_f1_thresh, test_aupr, test_cm, test_f1, test_mcc = evaluate_auc_gpu(
            Label, prob2)
        animator_1.add(epoch + 1, (None, test_loss.cpu().detach().numpy()), plotfile_loss)
        animator_2.add(epoch + 1, (None, test_auc), plotfile_auc)
        train_list.append(train_auc)
        test_list.append(test_auc)
        test_loss_list.append(test_loss)
        if len(test_loss_list) == 1:
            train_auc_1 = train_auc
            test_auc_1 = test_auc
            train_loss_1 = train_loss
            test_loss_1 = test_loss
            train_f1_ = train_f1
            test_f1_ = test_f1
            mcc = test_mcc
            cm = test_cm
            pro_ = prob2.cpu().detach().numpy()
            torch.save(net.state_dict(), mlp_model)
        else:
            if test_list[-1] > test_auc_1:
                train_auc_1 = train_auc
                test_auc_1 = test_auc
                train_loss_1 = train_loss
                test_loss_1 = test_loss
                train_f1_ = train_f1
                test_f1_ = test_f1
                mcc = test_mcc
                cm = test_cm
                pro_ = prob2.cpu().detach().numpy()
                torch.save(net.state_dict(), mlp_model)
    print(f'train loss {train_loss_1}, test loss {test_loss_1}, train auc '
          f'{train_auc_1:.3f}, test auc {test_auc_1:.3f}',
          f'train f1 {train_f1_}, test f1 {test_f1_}',
          f'mcc {mcc:.3f}, confusion_matrix {cm}')


def train_mse(net, train_iter, test_iter, loss, trainer, num_epochs, devices, 
              plotfile_loss, plotfile_auc, mlp_model):
    animator_1 = Animator(xlabel='epoch', xlim=[1, num_epochs],
                          legend=['train loss', 'test loss'])
    animator_2 = Animator(xlabel='epoch', xlim=[1, num_epochs],
                          legend=['train mae', 'test mae'])
    plotfile_mae = plotfile_auc
    test_list = []
    train_list = []
    test_loss_list = []
    for epoch in range(num_epochs):
        train_loss = 0
        for i, (features, labels) in enumerate(train_iter):
            l, pred = train_batch_ch13(net, features, labels, loss, trainer, devices)
            train_loss += l
            if i == 0:
                prob1 = pred
                Label = labels
            else:
                prob1 = torch.cat((prob1, pred), dim=0)
                Label = torch.cat((Label, labels), dim=0)
        train_loss = train_loss / (i+1)
        train_mae, train_r2 = evaluate_r2_gpu(Label, prob1)
        animator_1.add(epoch + 1, (train_loss.cpu().detach().numpy(), None), plotfile_loss)
        animator_2.add(epoch + 1, (train_mae, None), plotfile_mae)
        with torch.no_grad():
            net.eval()
            test_loss = 0
            for i, (features, labels) in enumerate(test_iter):
                X = features.to(devices[0])
                pred = net(X.float())
                y = labels.to(devices[0])
                test_loss += loss(torch.squeeze(pred, dim=1), y.float())
                if i == 0:
                    prob2 = pred
                    Label = labels
                else:
                    prob2 = torch.cat((prob2, pred), dim=0)
                    Label = torch.cat((Label, labels), dim=0)
            test_loss = test_loss / (i +1)
        test_mae, test_r2 = evaluate_r2_gpu(Label, prob2)
        animator_1.add(epoch + 1, (None, test_loss.cpu().detach().numpy()), plotfile_loss)
        animator_2.add(epoch + 1, (None, test_mae), plotfile_mae)
        train_list.append(train_mae)
        test_list.append(test_mae)
        if len(test_list) == 1:
            train_mae_1 = train_mae
            test_mae_1 = test_mae
            train_loss_1 = train_loss
            test_loss_1 = test_loss
            train_r2_1 = train_r2
            test_r2_1 = test_r2
            torch.save(net.state_dict(), mlp_model)
        else:
            if test_list[-1] < test_mae_1:
                train_mae_1 = train_mae
                test_mae_1 = test_mae
                train_loss_1 = train_loss
                test_loss_1 = test_loss
                train_r2_1 = train_r2
                test_r2_1 = test_r2
                torch.save(net.state_dict(), mlp_model)
    print(f'train loss {train_loss_1}, test loss {test_loss_1}, train mae '
          f'{train_mae_1:.3f}, test mae {test_mae_1:.3f}, train r2 '
          f'{train_r2_1:.3f}, test r2 {test_r2_1:.3f}')


def train_model(net, train_iter, test_iter, loss, trainer, num_epochs, devices,
                plotfile_loss, plotfile_auc, mlp_model):
    """Train a model with multiple GPUs.
    Defined in :numref:`sec_image_augmentation`"""
    num_batches = len(train_iter)
    # net = nn.DataParallel(net, device_ids=devices).to(devices[0])
    net = net.to(devices[0])
    if type(loss) is nn.MSELoss:
        train_mse(net, train_iter, test_iter, loss, trainer, num_epochs, devices, 
                  plotfile_loss, plotfile_auc, mlp_model)
    else:
        train_cls(net, train_iter, test_iter, loss, trainer, num_epochs, devices, 
                  plotfile_loss, plotfile_auc, mlp_model)


def init_weights(m):
    if isinstance(m, nn.Linear):
        nn.init.xavier_uniform_(m.weight)


# foca loss for banlance cross loss
class FocalLoss(nn.Module):

    def __init__(self, alpha=0.6, gamma=2, logits=True, reduction='mean', devices=None):
        super(FocalLoss, self).__init__()
        self.alpha = torch.tensor(alpha)
        self.gamma = torch.tensor(gamma)
        self.logits = logits
        self.reduction = reduction
        self.devices = devices

    def forward(self, inputs, targets):
        if self.logits:
            BCE_loss = F.binary_cross_entropy_with_logits(inputs,
                                                          targets.float(),
                                                          reduction='mean')
        else:
            BCE_loss = F.binary_cross_entropy(torch.sigmoid(inputs),
                                              targets.float(),
                                              reduction='mean')
        if inputs.is_cuda and not self.alpha.is_cuda:
            self.alpha = self.alpha.to(self.devices[0])
            self.gamma = self.gamma.to(self.devices[0])
        pt = torch.exp(-BCE_loss)
        F_loss = self.alpha * (1 - pt)**self.gamma * BCE_loss
        if self.reduction == 'mean':
            F_loss = F_loss.mean()
        elif self.reduction == 'sum':
            F_loss = F_loss.sum()
        return F_loss


def mlp(metadata,
        group,
        train_biom,
        test_biom,
        plotfile_loss,
        plotfile_auc,
        mlp_model,
        glove_embedding=None,
        loss=None,
        numb=1,
        p_drop=0.2,
        weight_decay=0.001,
        num_epochs=5,
        batch_size=128,
        num_hiddens=50,
        lr=0.0001):

    # with open(config_file) as f:
    #     config_defaults = yaml.load(f.read(),
    #                                 Loader=yaml.FullLoader)['parameters']
    # wandb.init(config=config_defaults, settings=wandb.Settings(start_method='fork'))
    # config = wandb.config
    # batch_size = config.batch_size
    # num_hiddens = config.num_hiddens
    # lr = config.lr

    if glove_embedding is not None:
        glove_embedding = f"{glove_embedding}_100.txt"
        feature_size = 100
    else:
        table = biom.load_table(train_biom)
        feature_size = table.shape[0]
    
    train_iter, test_iter, labels = load_data_imdb(batch_size, train_biom,
                                                   test_biom, metadata,
                                                   glove_embedding, group)

    devices = try_all_gpus(numb)
    net = Net(feature_size, num_hiddens, p_drop)
    
    # 训练
    trainer = torch.optim.Adam(net.parameters(), lr=lr, weight_decay=weight_decay)
    if loss is None:
        loss = nn.CrossEntropyLoss(reduction='mean')
    if loss == "FocalLoss":
        loss = FocalLoss(devices=devices)
    if loss == "BCE_loss":
        loss = nn.BCELoss(weight=None, reduction='mean')
    if loss == "MSE_loss":
        loss = nn.MSELoss(reduction='mean')
    train_model(net, train_iter, test_iter, loss, trainer, num_epochs, 
                devices, plotfile_loss, plotfile_auc, mlp_model)
