# 导入需要的库
import time
import torch
import biom
import numpy as np
import pandas as pd
from torch import nn

from matplotlib import pyplot as plt
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import confusion_matrix, f1_score, matthews_corrcoef, mean_absolute_error, r2_score, accuracy_score
from sklearn import metrics
from torch.nn import TransformerEncoder
import gc
import os
import matplotlib.pyplot as plt
import random


def gpu(i=0):
    """Select GPU device.
    """
    return torch.device(f'cuda:{i}')


def num_gpus():
    """Get the number of available GPUs.
    """
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    os.environ["CUDA_VISIBLE_DEVICES"] = "1"
    return torch.cuda.device_count()


def try_all_gpus(numb):
    """Return all available GPUs or [cpu()] if no GPU exists.
    """
    return [gpu(i) for i in [numb]]


def read_imdb(biom_table, metadata, labels_col, sample_id_col):
    """Read OTU table and sample metadata with flexible column mapping.
    
    Args:
        biom_table (str): Path to BIOM format OTU table
        metadata (str): Path to metadata file
        labels_col (str): Column name in metadata to use as labels (default: "group")
        sample_id_col (str, optional): Column name containing sample IDs. 
                                      Defaults to "sample_id"
    
    Returns:
        tuple: (normalized_otu_table, feature_ids, labels)
            - normalized_otu_table: Normalized OTU matrix (samples × features)
            - feature_ids: List of OTU feature IDs
            - labels: List of sample labels
    """
    labels = []
    table = biom.load_table(biom_table)
    feature_ids = table.ids(axis="observation")

    try:
        mapping_file = pd.read_csv(
            metadata,
            sep="\t",
            index_col=sample_id_col,  # 关键修改：使用用户指定的样本ID列
            dtype={sample_id_col: str},  # 强制样本ID为字符串类型
            low_memory=False
        )
    except ValueError as e:
        raise KeyError(f"Metadata missing required column: {sample_id_col}") from e

    ranked_table = table.rankdata(axis='sample', inplace=False)
    max_values = ranked_table.max(axis="sample")
    normalized_table = ranked_table.matrix_data.multiply(1 / max_values)
    normalized_table = normalized_table.toarray().T

    labels = mapping_file.loc[table.ids(axis='sample'), labels_col].tolist()

    return normalized_table, feature_ids, labels


class Fid:
    """Vocabulary mapping for OTU feature IDs with special token handling.

        Manages bidirectional mapping between OTU identifiers and integer indices,
        with built-in support for special tokens required by transformer models.
        Automatically handles unknown tokens and provides classification token.

    Attributes:
        idx_to_token (list): Ordered list where index corresponds to token
        token_to_idx (dict): Dictionary mapping tokens to indices
    
    Example:
        >>> otu_ids = ['OTU_123', 'OTU_456']
        >>> vocab = Fid(otu_ids, reserved_tokens=['<pad>'])
        >>> vocab['OTU_123']
        1  # Index order: ['cls',  'OTU_123', 'OTU_456', '<pad>', '<unk>',]
        >>> vocab.to_tokens(0)
        
    """

    def __init__(self, tokens, reserved_tokens=[]):

        # The list of tokens
        self.idx_to_token = list(
            sorted(
                set(reserved_tokens + ['<unk>'] + [
                    token
                    for token in tokens
                ])))
        # cls向量放置在第0号索引位置
        self.idx_to_token = ['cls'] + self.idx_to_token
        self.token_to_idx = {
            token: idx
            for idx, token in enumerate(self.idx_to_token)
        }

    def __len__(self):
        """Return number of tokens in vocabulary."""
        return len(self.idx_to_token)

    def __getitem__(self, tokens):
        if not isinstance(tokens, (list, tuple)):
            return self.token_to_idx.get(tokens, self.unk)
        return [self.__getitem__(token) for token in tokens]

    def to_tokens(self, indices):

        if hasattr(indices, '__len__') and len(indices) > 1:
            return [self.idx_to_token[int(index)] for index in indices]
        return self.idx_to_token[indices]

    @property
    def unk(self):
        return self.token_to_idx['<unk>']


class load_data_imdb(Dataset):
    """Dataset class for loading and processing OTU data."""

    def __init__(self, biom_table, metadata, labels_col, sample_id_col, num_steps=500):
        """Initialize OTU dataset processor.
        
        Args:
            biom_table (str): Path to BIOM-format table file (.biom)
            metadata (str): Path to sample metadata file (TSV format)
            labels_col (str): Metadata column name containing class labels
            sample_id_col (str): Metadata column containing sample IDs
            num_steps (int, optional): Maximum sequence length (excluding CLS token). 
                                     Default 500        """
        otu, fid, labels = read_imdb(biom_table, metadata, labels_col , sample_id_col)
        self.fid = Fid(fid, reserved_tokens=['<pad>'])
        self.features, self.abundance, self.mask = self.truncate_pad(otu, fid, num_steps)
        self.labels = torch.tensor(labels)
        self.otu = otu

    def truncate_pad(self, otu, fid, num_steps):
        """Truncate or pad OTU data to fixed length.
        
        Args:
            otu (ndarray): Raw OTU abundance matrix, shape [n_samples, n_features]
            fid (list): OTU feature IDs corresponding to columns in otu
            num_steps (int): Target sequence length (excluding CLS token)

        Returns:
            tuple: 
                - features (Tensor): Encoded sequences with CLS token, int64, [n_samples, num_steps+1]
                - abundance (Tensor): Abundance values, float32, [n_samples, num_steps+1]
                - mask (Tensor): Attention mask, int64, [n_samples, num_steps+1]
        """
        features = np.zeros((otu.shape[0], num_steps), dtype=int)
        cls_column = np.zeros((features.shape[0], 1)).astype(int)
        abundance = np.zeros((otu.shape[0], num_steps))
        cls_abundance = np.ones((features.shape[0], 1)).astype(int)
        for i in range(0, otu.shape[0]): 
            nonzero_count = np.count_nonzero(otu[i,])
            if nonzero_count >= num_steps:
                sorted_indices = np.argsort(otu[i,])[::-1]
                features[i, ] = np.array([self.fid[line]
                                          for line in
                                          fid[sorted_indices[:num_steps]]])
                abundance[i, ] = otu[i, sorted_indices[:num_steps]]
            else:
                nonzero_indices = np.nonzero(otu[i,])
                features[i, :nonzero_count] = np.array([self.fid[line]
                                                        for line in
                                                        fid[nonzero_indices]])
                features[i, nonzero_count:] = self.fid['<pad>']
                abundance[i, :nonzero_count] = otu[i, nonzero_indices]
        features = np.concatenate((cls_column, features), axis=1)
        abundance = np.concatenate((cls_abundance, abundance), axis=1)

        input_mask = torch.ones(features.shape[0], 
                                features.shape[1], 
                                dtype=torch.long)
        temp = features == self.fid['<pad>']
        input_mask[torch.tensor(temp)] = 0 # mask padding tokens
        input_mask[torch.tensor(features == 0)] = 0 # mask cls tokens
        return torch.tensor(features), torch.tensor(abundance), input_mask

    def __call__(self):
        """Get vocabulary object."""
        return self.fid

    def __getitem__(self, index):
        """Get sample by index.
        
        Args:
            index (int): Sample index.
            
        Returns:
            tuple: (features, abundance, label, mask)
        """
        return self.features[index, ], self.abundance[index, ], self.labels[index], self.mask[index]
    
    def __len__(self):
        """Get number of samples.
        
        Returns:
            int: Number of samples.
        """
        return self.otu.shape[0]


class TokenEmbedding:
    """Token embedding class for OTU features."""

    def __init__(self, data_dir):
        """Initialize embedding from file.
        
        Args:
            data_dir (str): Path to embedding file.
        """
        self.idx_to_token, self.idx_to_vec = self._load_embedding(data_dir)
        self.unknown_idx = 0
        self.token_to_idx = {
            token: idx
            for idx, token in enumerate(self.idx_to_token)
        }

    def _load_embedding(self, data_dir):
        """Load embedding vectors from file.
        
        Args:
            data_dir (str): Path to embedding file.
            
        Returns:
            tuple: (list of tokens, tensor of embedding vectors)
        """
        idx_to_token, idx_to_vec = ['<unk>'], []
        with open(data_dir, 'r') as f:
            for line in f:
                elems = line.rstrip().split(' ')
                token, elems = elems[0], [float(elem) for elem in elems[1:]]
                if len(elems) > 1:
                    idx_to_token.append(token)
                    idx_to_vec.append(elems)
        idx_to_vec = [[0] * len(idx_to_vec[0])] + idx_to_vec
        return idx_to_token, torch.tensor(idx_to_vec)

    def __getitem__(self, tokens):
        """Get embedding vectors for tokens.
        
        Args:
            tokens (list): List of tokens.
            
        Returns:
            torch.Tensor: Embedding vectors.
        """
        indices = [
            self.token_to_idx.get(token, self.unknown_idx) for token in tokens
        ]
        vecs = self.idx_to_vec[torch.tensor(indices)]
        return vecs

    def __len__(self):
        """Get number of tokens in vocabulary."""
        return len(self.idx_to_token)


class ScaledDotProductAttention(nn.Module):
    """Scaled Dot-Product Attention module."""

    def __init__(self, d_k):
        """Initialize attention module.
        
        Args:
            d_k (int): Dimension of key vectors.
        """
        super(ScaledDotProductAttention, self).__init__()
        self.d_k = d_k

    def forward(self, q, k, v, attn_mask):
        """Compute attention.
        
        Args:
            q (Tensor): Query tensor.
            k (Tensor): Key tensor.
            v (Tensor): Value tensor.
            attn_mask (Tensor): Attention mask tensor.
            
        Returns:
            tuple: (output tensor, attention weights)
        """
        attn_score = torch.matmul(q, k.transpose(-1, -2)) / np.sqrt(self.d_k)
        attn_score.masked_fill_(attn_mask, -1e9)
        attn_weights = nn.Softmax(dim=-1)(attn_score)
        output = torch.matmul(attn_weights, v)
        return output, attn_weights


class MultiHeadAttention(nn.Module):
    """Multi-head Attention module."""

    def __init__(self, d_model, n_heads):
        """Initialize multi-head attention.
        
        Args:
            d_model (int): Model dimension.
            n_heads (int): Number of attention heads.
        """
        super(MultiHeadAttention, self).__init__()
        self.n_heads = n_heads
        self.d_k = self.d_v = d_model // n_heads
        self.WQ = nn.Linear(d_model, d_model)
        self.WK = nn.Linear(d_model, d_model)
        self.WV = nn.Linear(d_model, d_model)
        self.scaled_dot_product_attn = ScaledDotProductAttention(self.d_k)
        self.linear = nn.Linear(n_heads * self.d_v, d_model)

    def forward(self, Q, K, V, attn_mask):
        """Compute multi-head attention.
        
        Args:
            Q (Tensor): Query tensor.
            K (Tensor): Key tensor.
            V (Tensor): Value tensor.
            attn_mask (Tensor): Attention mask tensor.
            
        Returns:
            tuple: (output tensor, attention weights)
        """
        batch_size = Q.size(0)
        q_heads = self.WQ(Q).view(batch_size, -1, self.n_heads,
                                  self.d_k).transpose(1, 2)
        k_heads = self.WK(K).view(batch_size, -1, self.n_heads,
                                  self.d_k).transpose(1, 2)
        v_heads = self.WV(V).view(batch_size, -1, self.n_heads,
                                  self.d_v).transpose(1, 2)
        attn_mask = attn_mask.unsqueeze(1).repeat(1, self.n_heads, 1, 1)
        attn, attn_weights = self.scaled_dot_product_attn(
            q_heads, k_heads, v_heads, attn_mask)
        attn = attn.transpose(1, 2).contiguous().view(batch_size, -1,
                                                      self.n_heads * self.d_v)
        output = self.linear(attn)
        return output, attn_weights


class PositionWiseFeedForwardNetwork(nn.Module):
    """Position-wise Feed Forward Network."""

    def __init__(self, d_model, d_ff):
        """Initialize FFN.
        
        Args:
            d_model (int): Model dimension.
            d_ff (int): Hidden layer dimension.
        """
        super(PositionWiseFeedForwardNetwork, self).__init__()
        self.linear1 = nn.Linear(d_model, d_ff)
        self.linear2 = nn.Linear(d_ff, d_model)
        self.relu = nn.ReLU()

    def forward(self, inputs):
        """Forward pass.
        
        Args:
            inputs (Tensor): Input tensor.
            
        Returns:
            Tensor: Output tensor.
        """
        output = self.relu(self.linear1(inputs))
        output = self.linear2(output)
        return output
    
    
class EncoderLayer(nn.Module):
    """Transformer Encoder Layer."""

    def __init__(self, d_model, n_heads, p_drop):
        """Initialize encoder layer.
        
        Args:
            d_model (int): Model dimension.
            n_heads (int): Number of attention heads.
            p_drop (float): Dropout probability.
        """
        super(EncoderLayer, self).__init__()
        self.mha = MultiHeadAttention(d_model, n_heads)
        self.dropout1 = nn.Dropout(p_drop)
        self.layernorm1 = nn.LayerNorm(d_model, eps=1e-6)
        self.ffn = PositionWiseFeedForwardNetwork(d_model, 4 * d_model)
        self.dropout2 = nn.Dropout(p_drop)
        self.layernorm2 = nn.LayerNorm(d_model, eps=1e-6)

    def forward(self, inputs, attn_mask):
        """Forward pass.
        
        Args:
            inputs (Tensor): Input tensor.
            attn_mask (Tensor): Attention mask tensor.
            
        Returns:
            tuple: (output tensor, attention weights)
        """
        attn_outputs, attn_weights = self.mha(inputs, inputs, inputs,
                                              attn_mask)
        attn_outputs = self.dropout1(attn_outputs)
        attn_outputs = self.layernorm1(inputs + attn_outputs)
        return attn_outputs, attn_weights


class TransformerEncoder(nn.Module):
    """Transformer Encoder for OTU data."""

    def __init__(self,
                 otu_size,
                 seq_len,
                 d_model=128,
                 n_layers=6,
                 n_heads=8,
                 p_drop=0.1,
                 d_ff=128,
                 pad_id=0):
        """Initialize transformer encoder.
        
        Args:
            otu_size (int): Size of OTU vocabulary.
            seq_len (int): Maximum sequence length.
            d_model (int, optional): Model dimension. Defaults to 128.
            n_layers (int, optional): Number of encoder layers. Defaults to 6.
            n_heads (int, optional): Number of attention heads. Defaults to 8.
            p_drop (float, optional): Dropout probability. Defaults to 0.1.
            d_ff (int, optional): Feed-forward dimension. Defaults to 128.
            pad_id (int, optional): Padding token ID. Defaults to 0.
        """
        super(TransformerEncoder, self).__init__()
        self.embedding = nn.Embedding(otu_size, d_model, padding_idx=pad_id)
        self.sinusoid_table = self.get_sinusoid_table(otu_size, d_model)
        self.pos_embedding = nn.Embedding.from_pretrained(self.sinusoid_table,
                                                          freeze=True)
        self.layers = nn.ModuleList([
            EncoderLayer(d_model, n_heads, p_drop)
            for _ in range(n_layers)
        ])
        self.sigmoid = nn.Sigmoid()
        self.relu = nn.ReLU()
        self.pad_id = pad_id
        self.dropout = nn.Dropout(p=p_drop)
        self.mlp = nn.Sequential(
                    nn.Linear(d_model, 1), 
                    )

    def forward(self, x, weight, mask, classification=True, encoder=False):
        """Forward pass.
        
        Args:
            x (Tensor): Input token IDs.
            weight (Tensor): Abundance weights.
            mask (Tensor): Attention mask.
            classification (bool, optional): Whether for classification task. Defaults to True.
            encoder (bool, optional): Whether to return encoder outputs. Defaults to False.
            
        Returns:
            tuple: (output tensor, attention weights) or (encoder outputs, embeddings)
        """
        inputs = self.embedding(x) + self.pos_embedding(x)
        inputs = inputs.permute(2, 0, 1) * weight
        inputs = inputs.permute(1, 2, 0)
        inputs_1 = inputs.clone()
        attn_pad_mask = self.get_attention_padding_mask(x, x, self.pad_id)

        for layer in self.layers:
            inputs, attn_weights = layer(inputs.float(), attn_pad_mask)
            
        embedding_sum  = inputs * mask.unsqueeze(-1)
        embedding = embedding_sum.sum(
            dim=1, keepdim=False) / mask.sum(1).unsqueeze(-1)
        outputs = self.dropout(embedding)
        outputs = self.relu(outputs)
        outputs = self.mlp(outputs.float())
        if classification is False:
            outputs = self.relu(outputs)
            
        if encoder == False:
            return outputs, attn_weights
        else:
            return inputs_1, embedding

    def get_attention_padding_mask(self, q, k, pad_id):
        """Create attention padding mask.
        
        Args:
            q (Tensor): Query tensor.
            k (Tensor): Key tensor.
            pad_id (int): Padding token ID.
            
        Returns:
            Tensor: Attention mask tensor.
        """
        attn_pad_mask = k.eq(pad_id).unsqueeze(1).repeat(1, q.size(1), 1)
        return attn_pad_mask

    def get_sinusoid_table(self, seq_len, d_model):
        """Create sinusoidal position encoding table.
        
        Args:
            seq_len (int): Sequence length.
            d_model (int): Model dimension.
            
        Returns:
            Tensor: Sinusoidal position encoding table.
        """
        def get_angle(pos, i, d_model):
            return pos / np.power(10000, (2 * (i // 2)) / d_model)

        sinusoid_table = np.zeros((seq_len, d_model))
        for pos in range(seq_len):
            for i in range(d_model):
                if i % 2 == 0:
                    sinusoid_table[pos, i] = np.sin(get_angle(pos, i, d_model))
                else:
                    sinusoid_table[pos, i] = np.cos(get_angle(pos, i, d_model))
        return torch.FloatTensor(sinusoid_table)



def evaluate_auc_gpu(Y, prob, device=None):
    """Evaluate AUC and other classification metrics.
    
    Args:
        Y (Tensor): True labels.
        prob (Tensor): Predicted probabilities.
        device (torch.device, optional): Device. Defaults to None.
        
    Returns:
        tuple: (AUC, optimal threshold, AUPR, confusion matrix, F1, MCC, accuracy)
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
    aupr = metrics.auc(recall, precision)
    y_hat = np.array([1 if i > max_f1_thresh else 0 for i in prob])
    cm = confusion_matrix(Y, y_hat)
    f1_scores = f1_score(Y, y_hat, average='macro')
    acc = accuracy_score(Y, y_hat)
    mcc = matthews_corrcoef(Y, y_hat)

    return auc, max_f1_thresh, aupr, cm, f1_scores, mcc, acc


def evaluate_r2_gpu(Y, pred, device=None):
    """Evaluate regression metrics.
    
    Args:
        Y (Tensor): True values.
        pred (Tensor): Predicted values.
        device (torch.device, optional): Device. Defaults to None.
        
    Returns:
        tuple: (MAE, R2 score)
    """
    Y = Y.numpy()
    pred = pred.cpu().detach().numpy()
    mae = mean_absolute_error(Y, pred)
    r2 = r2_score(Y, pred)

    return mae, r2


class Timer:
    """Timer for recording multiple running times."""

    def __init__(self):
        self.times = []
        self.start()

    def start(self):
        self.tik = time.time()

    def stop(self):
        self.times.append(time.time() - self.tik)
        return self.times[-1]

    def avg(self):
        return sum(self.times) / len(self.times)

    def sum(self):
        return sum(self.times)

    def cumsum(self):
        return np.array(self.times).cumsum().tolist()


def train_batch_ch13(net, X, y, abundance, loss, mask, trainer, devices):
    """Train for a minibatch with multiple GPUs.
    
    Args:
        net (nn.Module): 
            Neural network model. Should handle both classification and regression outputs.
            Expected to return tuple (predictions, attention_weights) from forward pass.
        X (Tensor or list): 
            Input features tensor(s). For single GPU: shape [batch_size, seq_len]
            If list, each element corresponds to sharded data for multi-GPU
        y (Tensor): 
            Target labels. Shape [batch_size] for classification, [batch_size, 1] for regression
        abundance (Tensor): 
            Relative abundance weights. Shape [batch_size, seq_len]
            Used for scaling embeddings in the model
        loss (nn.Module): 
            Loss criterion. Determines task type:
            - nn.MSELoss: Regression mode (direct output)
            - Other (e.g., BCELoss): Classification mode (sigmoid output)
        mask (Tensor): 
            Attention mask. Shape [batch_size, seq_len]
            1 = valid position, 0 = masked/padded
        trainer (optim.Optimizer): 
            Optimizer instance for parameter updates
        devices (list): 
            List of torch.device objects. Currently uses devices[0] for computation
        
    Returns:
        tuple: (loss sum, predictions)
    """
    if isinstance(X, list):
        X = [x.to(devices[0]) for x in X]
    else:
        X = X.to(devices[0])
    y = y.to(devices[0], dtype=torch.int64)
    abundance = abundance.to(devices[0])
    mask = mask.to(devices[0])
    net.train()
    trainer.zero_grad()
    if type(loss) is nn.MSELoss:
        pred, _ = net(X, abundance, classification=False)
        pred = torch.squeeze(pred, dim=1)
    else:
        pred, _ = net(X, abundance, mask)
        pred = torch.squeeze(pred, dim=1)
        pred = torch.sigmoid(pred)
    l = loss(pred.float(), y.float())
    l.backward()
    trainer.step()
    train_loss_sum = l.sum()
    return train_loss_sum, pred


def set_axes(axes, xlabel, ylabel, xlim, ylim, xscale, yscale, legend):
    """Set axes properties for matplotlib plots.
    
    Args:
        axes (matplotlib.axes.Axes): Axes object.
        xlabel (str): X-axis label.
        ylabel (str): Y-axis label.
        xlim (tuple): X-axis limits.
        ylim (tuple): Y-axis limits.
        xscale (str): X-axis scale.
        yscale (str): Y-axis scale.
        legend (list): Legend entries.
    """
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
                 figsize=(5, 5)):
        """Initialize animator.
        
        Args:
            xlabel (str, optional): X-axis label. Defaults to None.
            ylabel (str, optional): Y-axis label. Defaults to None.
            legend (list, optional): Legend entries. Defaults to None.
            xlim (tuple, optional): X-axis limits. Defaults to None.
            ylim (tuple, optional): Y-axis limits. Defaults to None.
            xscale (str, optional): X-axis scale. Defaults to 'linear'.
            yscale (str, optional): Y-axis scale. Defaults to 'linear'.
            fmts (tuple, optional): Line formats. Defaults to ('-', 'm--', 'g-.', 'r:').
            nrows (int, optional): Number of rows. Defaults to 1.
            ncols (int, optional): Number of columns. Defaults to 1.
            figsize (tuple, optional): Figure size. Defaults to (5, 5).
        """
        if legend is None:
            legend = []
        self.fig, self.axes = plt.subplots(nrows, ncols, figsize=figsize)
        if nrows * ncols == 1:
            self.axes = [
                self.axes,
            ]
        self.config_axes = lambda: set_axes(self.axes[0], xlabel, ylabel, xlim,
                                            ylim, xscale, yscale, legend)
        self.X, self.Y, self.fmts = None, None, fmts

    def add(self, x, y, plotfile):
        """Add data points to the figure.
        
        Args:
            x: X values.
            y: Y values.
            plotfile (str): Path to save plot.
        """
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
        
        
def train_model(net, train_iter, test_iter, loss,trainer, scheduler, num_epochs, devices,
                embedding_birnn, plotfile_loss, plotfile_auc):
    """Train model with multiple GPUs.
    
    Args:
        net (nn.Module): Neural network model.
        train_iter (DataLoader): Training data loader.
        test_iter (DataLoader): Test data loader.
        loss (nn.Module): Loss function.
        trainer (optim.Optimizer): Optimizer.
        scheduler (optim.lr_scheduler): Learning rate scheduler.
        num_epochs (int): Number of training epochs.
        devices (list): List of devices.
        embedding_birnn (str): Path to save model.
        plotfile_loss (str): Path to save loss plot.
        plotfile_auc (str): Path to save AUC plot.
    """
    if type(loss) is nn.MSELoss:
        train_regression(net, train_iter, test_iter, loss,trainer, scheduler, num_epochs, devices,
                embedding_birnn, plotfile_loss, plotfile_auc)
    else:
        train_cls(net, train_iter, test_iter, loss,trainer, scheduler, num_epochs, devices,
                embedding_birnn, plotfile_loss, plotfile_auc)


def train_cls(net, train_iter, test_iter, loss,trainer, scheduler, num_epochs, devices,
                embedding_birnn, plotfile_loss, plotfile_auc):
    """Transformer-based Classification Model for Microbiome Data Analysis
    
    Args:
        net (nn.Module): 
            Neural network model to be trained. Should output logits for classification.
        train_iter (DataLoader): 
            Training data loader yielding tuples of (features, abundance, labels, mask).
            - features: Tensor of shape [batch_size, seq_len] containing OTU indices
            - abundance: Tensor of shape [batch_size, seq_len] with relative abundances
            - labels: Tensor of shape [batch_size] with class labels (0/1)
            - mask: Tensor of shape [batch_size, seq_len] for attention masking
        test_iter (DataLoader): 
            Validation/Test data loader with same structure as train_iter
        loss (nn.Module): 
            Loss function module (e.g., nn.BCELoss). Should handle sigmoid outputs
        trainer (optim.Optimizer): 
            Optimizer instance (e.g., Adam)
        scheduler (lr_scheduler._LRScheduler): 
            Learning rate scheduler (e.g., ExponentialLR)
        num_epochs (int): 
            Maximum number of training epochs
        devices (list): 
            List of torch.device objects for GPU acceleration (single GPU supported)
        embedding_birnn (str): 
            Path to save best model weights (state_dict format)
        plotfile_loss (str): 
            Path to save loss curve plot (PNG format)
        plotfile_auc (str): 
            Path to save AUC curve plot (PNG format)
    Returns:
        None: 
            Model weights are saved to disk. Training metrics printed to stdout.
    """
    plotfile_auc = plotfile_auc
    animator_1 = Animator(xlabel='epoch', xlim=[1, num_epochs],legend=['train loss', 'test loss'])
    animator_2 = Animator(xlabel='epoch', xlim=[1, num_epochs], ylim=[0, 1],
                        legend=['train auc', 'test auc'])
    net = net.to(devices[0])
    test_list = []
    train_list = []
    test_loss_list = []
    
    for epoch in range(num_epochs):
        train_loss = 0
        for i, (features, abundance, labels, mask) in enumerate(train_iter):
            l, pred = train_batch_ch13(
                net, features, labels, abundance, loss, mask, trainer, devices)
            train_loss += l
            if i == 0:
                prob1 = pred
                Label = labels
            else:
                prob1 = torch.cat((prob1, pred), dim=0)
                Label = torch.cat((Label, labels), dim=0)
                
        train_loss /= (i + 1)
        animator_1.add(epoch + 1, (train_loss.cpu().detach().numpy(), None), plotfile_loss)
        train_auc, train_f1_thresh, train_aupr, train_cm, train_f1, train_mcc, train_acc = evaluate_auc_gpu(
            Label, prob1)
        animator_2.add(epoch + 1, (train_auc, None), plotfile_auc)
        
        with torch.no_grad():
            net.eval()
            test_loss = 0
            for i, (features, abundance, labels, mask) in enumerate(test_iter):
                X = features.to(devices[0])
                abundance = abundance.to(devices[0])
                mask = mask.to(devices[0])
                pred, _ = net(X, abundance, mask)
                pred = torch.squeeze(pred, dim=1)
                y = labels.to(devices[0], dtype=torch.int64)
                pred = torch.sigmoid(pred)
                l = loss(pred.float(), y.float())
                test_loss += l
                if i == 0:
                    prob2 = pred
                    Label = labels
                else:
                    prob2 = torch.cat((prob2, pred), dim=0)
                    Label = torch.cat((Label, labels), dim=0)
            test_loss /= (i + 1)
            
        animator_1.add(epoch + 1, (None, test_loss.cpu().detach().numpy()), plotfile_loss)
        test_auc, test_f1_thresh, test_aupr, test_cm, test_f1, test_mcc, test_acc = evaluate_auc_gpu(
            Label, prob2)
        animator_2.add(epoch + 1, (None, test_auc), plotfile_auc)
        train_list.append(train_auc)
        test_list.append(test_auc)
        test_loss_list.append(test_loss)
        
        if len(test_loss_list) == 1:
            train_auc_1 = train_auc
            test_auc_1 = test_auc
            train_loss_1 = train_loss
            test_loss_1 = test_loss
            train_aupr_1 = train_aupr
            test_aupr_1 = test_aupr
            train_f1_1 = train_f1
            test_f1_1 = test_f1
            train_mcc_1 = train_mcc
            test_mcc_1 = test_mcc
            train_cm_1 = train_cm
            test_cm_1 = test_cm
            train_acc_1 = train_acc
            test_acc_1 = test_acc
            pro_ = prob2.cpu().detach().numpy()
            torch.save(net.state_dict(), embedding_birnn)
        else:
            if test_list[-1] > test_auc_1:
                train_auc_1 = train_auc
                test_auc_1 = test_auc
                train_loss_1 = train_loss
                test_loss_1 = test_loss
                train_aupr_1 = train_aupr
                test_aupr_1 = test_aupr
                train_f1_1 = train_f1
                test_f1_1 = test_f1
                train_mcc_1 = train_mcc
                test_mcc_1 = test_mcc
                train_cm_1 = train_cm
                test_cm_1 = test_cm
                train_acc_1 = train_acc
                test_acc_1 = test_acc
                pro_ = prob2.cpu().detach().numpy()
                torch.save(net.state_dict(), embedding_birnn)
                
    print(f'train loss {train_loss_1}, test loss {test_loss_1}, train auc',
        f'{train_auc_1:.3f}, test auc {test_auc_1:.3f}',
        f'train f1 {train_f1_1:.3f}, test f1 {test_f1_1:.3f}',
        f'train mcc {train_mcc_1:.3f}, test mcc {test_mcc_1:.3f}',
        f'train aupr {train_aupr_1:.3f}, test aupr {test_aupr_1:.3f}',
        f'train acc {train_acc_1:.3f}, test acc {test_acc_1:.3f}',
        f'confusion_matrix {test_cm_1}')
    del net, features, labels, prob1, prob2, Label
    gc.collect()
    torch.cuda.empty_cache()


def train_regression(net, train_iter, test_iter, loss,trainer, scheduler, num_epochs, devices,
                embedding_birnn, plotfile_loss, plotfile_auc):
    """Train regression model.
    
    Args:
        net (nn.Module): Neural network model.
        train_iter (DataLoader): Training data loader.
        test_iter (DataLoader): Test data loader.
        loss (nn.Module): Loss function.
        trainer (optim.Optimizer): Optimizer.
        scheduler (optim.lr_scheduler): Learning rate scheduler.
        num_epochs (int): Number of training epochs.
        devices (list): List of devices.
        embedding_birnn (str): Path to save model.
        plotfile_loss (str): Path to save loss plot.
        plotfile_auc (str): Path to save MAE plot.
    """
    plotfile_mae = plotfile_auc
    animator_1 = Animator(xlabel='epoch', xlim=[1, num_epochs],legend=['train loss', 'test loss'])
    animator_2 = Animator(xlabel='epoch', xlim=[1, num_epochs],legend=['train mae', 'test mae'])
    net = net.to(devices[0])
    test_list = []
    train_list = []
    test_loss_list = []
    true_age = []
    pred_age = []
    
    for epoch in range(num_epochs):
        train_loss = 0
        for i, (features, abundance, labels, mask) in enumerate(train_iter):
            l, pred = train_batch_ch13(
                net, features, labels, abundance, loss, mask,trainer, devices)
            train_loss = train_loss + l
            if i == 0:
                prob1 = pred
                Label = labels
            else:
                prob1 = torch.cat((prob1, pred), dim=0)
                Label = torch.cat((Label, labels), dim=0)
                
        train_loss = train_loss / (i + 1)
        animator_1.add(epoch + 1, (train_loss.cpu().detach().numpy(), None), plotfile_loss)
        train_mae, train_r2= evaluate_r2_gpu(Label, prob1)
        animator_2.add(epoch + 1, (train_mae, None), plotfile_mae)
        
        true_age_test=[]
        pred_age_test=[]
        with torch.no_grad():
            net.eval()
            test_loss = 0
            for i, (features, abundance, labels, mask) in enumerate(test_iter):
                X = features.to(devices[0])
                abundance = abundance.to(devices[0])
                mask = mask.to(devices[0])
                pred, _ = net(X, abundance, mask, classification=False)
                pred = torch.squeeze(pred, dim=1)
                y = labels.to(devices[0], dtype=torch.int64)
                
                # true_age_test.append(labels.detach().cpu().numpy())
                # pred_age_test.append(pred.detach().cpu().numpy())
                
                l = loss(pred.float(), y.float())
                test_loss += l
                if i == 0:
                    prob2 = pred
                    Label = labels
                else:
                    prob2 = torch.cat((prob2, pred), dim=0)
                    Label = torch.cat((Label, labels), dim=0)
                    
            # true_age_test=np.concatenate(true_age_test)
            # pred_age_test=np.concatenate(pred_age_test)
            # true_age.append(true_age_test)
            # pred_age.append(pred_age_test)
            test_loss /= (i + 1)
            
        animator_1.add(epoch + 1, (None, test_loss.cpu().detach().numpy()), plotfile_loss)
        test_mae, test_r2 = evaluate_r2_gpu(Label, prob2)
        animator_2.add(epoch + 1, (None, test_mae), plotfile_mae)
        train_list.append(train_mae)
        test_list.append(test_mae)
        test_loss_list.append(test_loss)
        
        if len(test_loss_list) == 1:
            train_r2_1 = train_r2
            test_r2_1 = test_r2
            train_loss_1 = train_loss
            test_loss_1 = test_loss
            train_mae_1 = train_mae
            test_mae_1 = test_mae
            pro_ = prob2.cpu().detach().numpy()
            torch.save(net.state_dict(), embedding_birnn)
        else:
            if test_list[-1] < test_mae_1:
                train_r2_1 = train_r2
                test_r2_1 = test_r2
                train_loss_1 = train_loss
                test_loss_1 = test_loss
                train_mae_1 = train_mae
                test_mae_1 = test_mae
                torch.save(net.state_dict(), embedding_birnn)
                
    print(f'train loss {train_loss_1}, test loss {test_loss_1}, train mae '
        f'{train_mae_1:.3f}, test mae {test_mae_1:.3f}, train r2 '
        f'{train_r2_1:.3f}, test r2 {test_r2_1:.3f}')
    del net, features, labels, prob1, prob2, Label
    gc.collect()
    torch.cuda.empty_cache()


def init_transformer_weights(model, method='xavier', exclude='embedding'):
    """Initialize transformer weights.
    
    Args:
        model (nn.Module): Model to initialize.
        method (str, optional): Initialization method. Defaults to 'xavier'.
        exclude (str, optional): Parameter names to exclude. Defaults to 'embedding'.
    """
    for name, w in model.named_parameters():
        if exclude not in name:
            if 'weight' in name:
                if method == 'xavier':
                    if len(w.shape) < 2:
                        nn.init.xavier_normal_(w.unsqueeze(0))
                    else:
                        nn.init.xavier_normal_(w)
                elif method == 'kaiming':
                    if len(w.shape) < 2:
                        nn.init.kaiming_normal_(w.unsqueeze(0))
                    else:
                        nn.init.kaiming_normal_(w)
                else:
                    nn.init.normal_(w)
            elif 'bias' in name:
                nn.init.constant_(w, 0)
            else:
                pass


class FocalLoss(nn.Module):
    """Focal loss for imbalanced classification."""

    def __init__(self, alpha=0.6, gamma=2, reduction='mean', devices=None):
        """Initialize focal loss.
        
        Args:
            alpha (float, optional): Weighting factor. Defaults to 0.6.
            gamma (float, optional): Focusing parameter. Defaults to 2.
            reduction (str, optional): Reduction method. Defaults to 'mean'.
            devices (list, optional): List of devices. Defaults to None.
        """
        super(FocalLoss, self).__init__()
        self.alpha = torch.tensor(alpha)
        self.gamma = torch.tensor(gamma)
        self.reduction = reduction
        self.devices = devices

    def forward(self, inputs, targets):
        """Compute focal loss.
        
        Args:
            inputs (Tensor): Predicted probabilities.
            targets (Tensor): Target labels.
            
        Returns:
            Tensor: Computed loss.
        """
        if inputs.is_cuda and not self.alpha.is_cuda:
            self.alpha = self.alpha.to(self.devices[0])
            self.gamma = self.gamma.to(self.devices[0])
        pt = inputs
        alpha = self.alpha
        F_loss = - alpha * (1 - pt) ** self.gamma * targets * torch.log(pt) - \
               (1 - alpha) * pt ** self.gamma * (1 - targets) * torch.log(1 - pt)
        if self.reduction == 'mean':
            F_loss = torch.mean(F_loss)
        elif self.reduction == 'sum':
            F_loss = torch.sum(F_loss)
        return F_loss

def set_seed(seed=11):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    

def Attention_biom(
        metadata,
        train_biom,
        test_biom,
        embedding_birnn,
        plotfile_loss,
        plotfile_auc,
        labels_col="group",
        sample_id_col="sample_id",
        num_steps=400,
        p_drop=0,
        d_ff=64,
        batch_size=128,
        d_model=100,
        n_layers=2,
        n_heads=2,
        numb=1,
        lr=0.0005,
        weight_decay=0,
        num_epochs=100,
        loss="BCE_loss",
        alpha=0.6,
        glove_embedding=None,):
    """End-to-end training pipeline for OTU-based microbial analysis using transformer models.
    
    Args:
      metadata (str): 
            Path to TSV metadata file containing sample labels. Expected format:
            - Must contain columns: <sample_id_col>, <labels_col>
            - Sample IDs must match BIOM table entries
        train_biom (str): 
            Path to training BIOM format file (features x samples)
        test_biom (str): 
            Path to test/validation BIOM file (same format as train_biom)
        embedding_birnn (str): 
            Output path for saving best model weights (.pt format)
        plotfile_loss (str): 
            Output path for loss curve plot (e.g., 'loss.png')
        plotfile_auc (str): 
            Output path for AUC/MAE curve plot (e.g., 'auc.png')
        labels_col (str, optional): 
            Metadata column name containing classification labels. Default "group"
        sample_id_col (str, optional): 
            Metadata column with sample IDs matching BIOM. Default "sample_id"
        num_steps (int, optional): 
            Maximum OTU sequence length (after CLS token). Default 400
        p_drop (float, optional): 
            Dropout probability for transformer layers. Default 0 (no dropout)
        d_ff (int, optional): 
            Feed-forward dimension in transformer. Default 64
        batch_size (int, optional): 
            Training batch size. Default 128
        d_model (int, optional): 
            Transformer embedding dimension. Default 100
        n_layers (int, optional): 
            Number of transformer encoder layers. Default 2
        n_heads (int, optional): 
            Number of attention heads per layer. Default 2
        numb (int, optional): 
            GPU device index (0-based). Default 1 (cuda:1)
        lr (float, optional): 
            Learning rate for Adam optimizer. Default 0.0005
        weight_decay (float, optional): 
            L2 regularization strength. Default 0 (no regularization)
        num_epochs (int, optional): 
            Maximum training epochs. Default 100
        loss (str, optional): 
            Loss function specification. Options:
            - "BCE_loss": Binary Cross Entropy (classification)
            - "FocalLoss": Focal Loss for class imbalance
            - "MSE_loss": Mean Squared Error (regression)
            - "MAE_loss": Mean Absolute Error (regression)
            Default "BCE_loss"
        alpha (float, optional): 
            Alpha parameter for Focal Loss. Default 0.6
        glove_embedding (str, optional): 
            Path to GloVe embeddings file for OTU initialization. 
    """
    set_seed(11)
    devices = try_all_gpus(numb)

    # if glove_embedding is not None:
    #     glove_embedding = f"{glove_embedding}_{d_model}.txt"

    train_data = load_data_imdb(train_biom,
                                metadata,
                                labels_col,
                                sample_id_col,
                                num_steps)
    test_data = load_data_imdb(test_biom,
                               metadata,
                               labels_col,
                               sample_id_col,
                               num_steps)
    train_iter = DataLoader(train_data,
                            batch_size=batch_size,
                            shuffle=True)
    test_iter = DataLoader(test_data,
                           batch_size=batch_size,
                           shuffle=False)

    fid_dict = train_data()
    net = TransformerEncoder(otu_size=len(fid_dict),
                            seq_len=num_steps+1,
                            d_model=d_model,
                            n_layers=n_layers,
                            n_heads=n_heads,
                            p_drop=p_drop,
                            d_ff=d_ff,
                            pad_id=fid_dict['<pad>'])
    net.apply(init_transformer_weights)

    if glove_embedding is not None:
        glove_embedding = TokenEmbedding(glove_embedding)
        embeds = glove_embedding[fid_dict.idx_to_token]
        net.embedding.weight.data.copy_(embeds)
    else:
        net.embedding.weight.data = torch.zeros(len(fid_dict), 100)
    net.embedding.weight.requires_grad = False

    trainer = torch.optim.Adam(
        net.parameters(), lr=lr, weight_decay=weight_decay)
    scheduler = torch.optim.lr_scheduler.ExponentialLR(trainer, gamma=0.9)
    if loss is None:
        loss = nn.CrossEntropyLoss(reduction='mean')
    if loss == "FocalLoss":
        loss = FocalLoss(alpha=alpha, devices=devices)
    if loss == "BCE_loss":
        loss = nn.BCELoss(weight=None, reduction='mean')
    if loss == "MAE_loss":
        loss = nn.L1Loss(reduction='mean')
    if loss == "MSE_loss":
        loss = nn.MSELoss(reduction='mean')
        
    train_model(net, train_iter, test_iter, loss, trainer, scheduler, num_epochs, devices,
                embedding_birnn, plotfile_loss, plotfile_auc)