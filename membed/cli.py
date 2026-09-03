# -*- coding: utf-8 -*-
"""Console script for membed."""
import sys
from functools import partial
import click

from .cooccur_embedding import (SUPPORTED_METRICS, get_feature_dict,
                                cooccur_workflow, train_glove_model,
                                build_x_max_file_workflow)
from .util import set_log_level
from .otu_attention import Attention_biom


click.option = partial(click.option, show_default=True)


def decorator_composer(*decorators):

    def deco(f):
        for dec in reversed(decorators):
            f = dec(f)
        return f

    return deco

_common_parameters = decorator_composer(
    click.option(
        '-f',
        '--format',
        type=click.Choice(['biom', 'tsv']),
        default='biom',
        help=
        'Input table format. For the tsv table, the features are in row and samples are in column; '
        'columns are separated by TABs; the 1st column must be feature IDs and the 1st row must be sample IDs; '
        'lines starting with `#` will be ignored.'),
    click.option('-l',
                 '--log',
                 type=click.Choice([
                     'critical', 'error', 'warning', 'info', 'debug', 'notset'
                 ]),
                 default='info',
                 help='Set logging level.'),
    click.option('--force',
                 is_flag=True,
                 help='Overwrite output file if it already exists.'))


@click.group(context_settings=dict(
    help_option_names=['-h', '--help'],
    # allow case insensitivity for the (sub)commands and their options
    token_normalize_func=lambda x: x.lower()))
def main():
    """Insight into microbial embedding."""


@main.command()
@click.option('-b',
              '--biom-file',
              type=click.Path(exists=True),
              required=True,
              help='Input BIOM-format feature table file path')
@click.option('-d',
              '--feature-dict',
              type=click.STRING,
              required=True,
              help='Output feature dictionary file path. Format: two columns [feature_id] [non_zero_count], space-separated.')
@_common_parameters
def dict(**kwargs):
    """Count the number of non-zero occurrences for each feature across all samples in the dataset.
    Parameters:
        biom_file (str): Path to input BIOM format file
        feature_dict (str): Output path for feature-count dictionary
        log (str): Logging level (default: info)
        force (bool): Overwrite existing files

    Example:
        $ membed dict -b table.biom -d feature-dict.csv
    """
    set_log_level(kwargs['log'])
    params = {key: kwargs[key] for key in ('biom_file', 'feature_dict')}
    get_feature_dict(**params)


@main.command()
@click.option('-b',
              '--biom-file',
              type=click.Path(exists=True),
              required=True,
              help='Input BIOM-format feature table path.')
@click.option('-c',
              '--cooccur-file',
              type=click.STRING,
              required=True,
              help='Output co-occurrence record file path. Headerless binary array of '
              'dtype [(otu1, i4), (otu2, i4), (value, f8)], symmetrically expanded, '
              'with 1-based feature indices (GloVe reserves index 0).')
@click.option('--metric',
              type=click.Choice(SUPPORTED_METRICS),
              default='russell_rao',
              help='''Similarity metric (default: russell_rao).
All metrics are similarities: higher = more co-occurring.

1. **russell_rao** - Binary presence/absence similarity
   - Formula: (A∩B) / TotalSamples
   - Use Case: Basic co-occurrence frequency
   - Requires: Binary data (auto-converted)
   - Example: `--metric russell_rao`

2. **weighted_russell_rao** - Rank-distance weighted co-occurrence
   - Formula: ∑ 1/|rank_u - rank_v| over co-occurring samples / TotalSamples
   - Use Case: Embedding-focused relationships (Eq. 2)
   - Requires: Within-sample ranks (auto-computed)
   - Example: `--metric weighted_russell_rao --cpus 20`

3. **jaccard** - Jaccard similarity coefficient
   - Formula: (A∩B)/(A∪B)
   - Range: 0-1
   - Requires: Binary data (auto-converted)
   - Example: `--metric jaccard`

4. **faith** - Faith's similarity
   - Formula: (A∩B + (¬A∩¬B)/n)/n
   - Requires: Binary data (auto-converted)
   - Example: `--metric faith`

5. **braycurtis_percentile** / **braycurtis_totalsum** - Bray-Curtis similarity
   - Formula: ∑min(u,v)/(∑u+∑v); dissimilarity is 1 - 2x this
   - Suffix picks the scaling: percentile = within-sample rank / max,
     totalsum = relative abundance
   - Example: `--metric braycurtis_totalsum`

6. **abundance_percentile** / **abundance_totalsum** - Abundance-weighted co-occurrence
   - Formula: mean(min(u,v) * (1 - |u-v|))
   - Suffix picks the scaling, as above
   - Example: `--metric abundance_percentile --cpus 20`'''
)
@click.option('--cpus',
              type=click.INT,
              default=1,
              help='Number of CPU cores for the pairwise scan. Ignored by '
              'russell_rao, which is a single matrix product.')
@_common_parameters
def cooccur(**kwargs):
    """Compute pairwise cooccurrence of microbes in a biom table.

    Example:
    $ membed cooccur -b table.biom -c table.co  --metric russell_rao
    $ membed cooccur -b table.biom -c table.co  --metric abundance_percentile --cpus 20
    """
    set_log_level(kwargs['log'])
    params = {
        key: kwargs[key]
        for key in ('biom_file', 'cooccur_file',  'metric', 'cpus' )
    }
    cooccur_workflow(**params)


@main.command()
@click.option('-c',
              '--cooccur-file',
              type=click.STRING,
              required=True,
              help='Input co-occurrence record file written by `membed cooccur`.')
@click.option('-x',
              '--x-max-file',
              type=click.STRING,
              required=True,
              help="Output x_max value file path; '.npy' is appended when missing. "
              "Used to truncate the GloVe weighting function.")
@click.option('--percentile-num',
              '--percentile_num',
              'percentile_num',
              type=click.FloatRange(0, 100),
              default=80,
              help='Percentile (0-100) of the co-occurrence values to use as x_max. '
              'This study adopts the 80th percentile.')
@_common_parameters
def build_x_max_file(**kwargs):
    """Pick the x_max weighting cutoff from the co-occurrence records.

    Example:
    $ membed build-x-max-file -c table.co -x xmax_file.npy --percentile-num 80
    """
    set_log_level(kwargs['log'])

    params = {
        key: kwargs[key]
        for key in ('cooccur_file', 'x_max_file', 'percentile_num' )
    }

    build_x_max_file_workflow(**params)


@main.command()
@click.option('-d',
              '--feature-dict',
              type=click.STRING,
              required=True,
              help='Input feature dictionary written by `membed dict`, used as '
              "GloVe's vocab file.")
@click.option('-c',
              '--cooccur-file',
              type=click.STRING,
              required=True,
              help='Input co-occurrence record file written by `membed cooccur`.')
@click.option('-x',
              '--x-max-file',
              type=click.STRING,
              required=True,
              help='Input x_max file written by `membed build-x-max-file`; a missing '
              "'.npy' suffix is added automatically.")
@click.option('-r',
              '--result',
              type=click.STRING,
              required=True,
              help='Output directory, created if absent. Holds '
              'embeddings_<embedding-size>.txt plus the gradsq and shuffle files.')
@click.option('--lr',
              type=click.FLOAT,
              default=0.05,
              help='Learning rate (eta) for AdaGrad.')
@click.option('--embedding-size',
              type=click.INT,
              default=100,
              help='Dimension of the output embeddings.')
@click.option('--iter',
              type=click.INT,
              default=50,
              help='Number of training iterations over the co-occurrence records.')
@click.option('--cpus',
              type=click.INT,
              default=4,
              help='Number of threads for the GloVe binary.')
@_common_parameters
def glove_train(**kwargs):
    """Run GloVe embedding.

    Example:
    $ membed glove-train -d feature-dict.csv -c table.co -r ./result/ -x xmax_file.npy --lr 0.05 --embedding-size 100 --iter 50 --cpus 2
    """
    set_log_level(kwargs['log'])
    params = {
        key: kwargs[key]
        for key in ("cooccur_file", "feature_dict", "result",
                     "x_max_file", "lr", "embedding_size",
                    "cpus")
    }
    # The CLI flag stays --iter; the function avoids shadowing the builtin.
    params["n_iter"] = kwargs["iter"]
    train_glove_model(**params)

@main.command()
@click.option('-tra_otu',
              '--train-biom',
              type=click.STRING,
              required=True,
              help='Path to training dataset in BIOM format')
@click.option('-val_otu',
              '--valid-biom',
              type=click.STRING,
              required=True,
              help='Path to validation dataset in BIOM format. Drives early '
              'stopping, checkpoint selection and the decision threshold; never '
              'trained on.')
@click.option('-tes_otu',
              '--test-biom',
              type=click.STRING,
              required=True,
              help='Path to testing dataset in BIOM format')
@click.option('-m',
              '--metadata',
              type=click.STRING,
              required=True,
              help='Path to a TAB-separated metadata file covering the samples of '
              'both tables.')
@click.option('--labels-col',
              '--labels_col',
              'labels_col',
              type=click.STRING,
              required=True,
              help='Metadata column holding the class labels. Must already be '
              'numeric 0/1; encode the classes before training.')
@click.option('--sample-id-col',
              '--sample_id_col',
              'sample_id_col',
              type=click.STRING,
              required=False,
              default='sample_id',
              help="Unique sample ID column name matching BIOM sample IDs")
@click.option('-e',
              '--embedding-birnn',
              type=click.STRING,
              required=True,
              help="Output path for the best epoch's model state_dict.")
@click.option('-ploss',
              '--plotfile-loss',
              type=click.STRING,
              required=True,
              help='Output path for training loss plot')
@click.option('-pauc',
              '--plotfile-auc',
              type=click.STRING,
              required=True,
              help='Output path for the train/validation ROC-AUC curve plot')
@click.option('--num-steps',
              type=click.INT,
              default=400,
              help='Number of OTU positions kept per sample (truncated or padded).')
@click.option('--p-drop',
              type=click.FLOAT,
              default=0.0,
              help='Dropout probability')
@click.option('--batch-size',
              type=click.INT, 
              default=64,
              help='Training batch size')
@click.option('--d-model', 
              type=click.INT, 
              default=100,
              help='Model embedding dimension')
@click.option('--d-ff',
              type=click.INT,
              default=None,
              help='Width of the feed-forward hidden layer in each encoder layer. '
              'Defaults to 4 * --d-model.')
@click.option('--head-hidden',
              type=click.INT,
              default=None,
              help='Width of the hidden layer in the output head. Defaults to '
              'd_model // 2; pass 0 for the single-linear-layer head.')
@click.option('--n-layers', 
              type=click.INT, 
              default=2,
              help='Number of transformer layers')
@click.option('--n-heads', 
              type=click.INT, 
              default=2,
              help='Number of attention heads')
@click.option('--numb',
              type=click.INT,
              default=1,
              help='Zero-based CUDA device index to train on; negative values '
              'run on CPU.')
@click.option('--lr', 
              type=click.FLOAT, 
              default=0.0005,
              help='Learning rate')
@click.option('--weight-decay', 
              type=click.FLOAT, 
              default=0.0,
              help='L2 regularization weight')
@click.option('--num-epochs', 
              type=click.INT, 
              default=1,
              help='Number of training epochs')
@click.option('--loss',
              type=click.Choice(['BCE_loss', 'BCEWithLogits', 'FocalLoss']),
              default='BCE_loss',
              help='Criterion to train with.')
@click.option('--alpha',
              type=click.FLOAT,
              default=0.6,
              help='Positive-class weight, used only by FocalLoss.')
@click.option('-g',
              '--glove-embedding',
              type=click.STRING,
              default=None,
              help='Path to the pretrained embedding file written by '
              '`membed glove-train` (e.g. result/embeddings_100.txt). Its vector '
              'width must equal --d-model. Omitted, the frozen embedding table is '
              'filled with fixed random codes instead.')
@click.option('--pred-out',
              type=click.STRING,
              default=None,
              help='Prefix for the prediction files, written as <pred-out>_valid.csv '
              'and <pred-out>_test.csv. Defaults to --embedding-birnn without its '
              'extension.')
@_common_parameters
def class_attention(**kwargs):
    """Run the attention classifier on a train/test pair of BIOM tables.

    Training, validation and test come in as three separate tables; the epoch to
    keep and the decision threshold are chosen on validation only. The test
    table is scored once at the end with the frozen threshold, and per-sample
    predictions are written next to the model.

    Example:
        membed class-attention \\
            --train-biom train.biom \\
            --valid-biom valid.biom \\
            --test-biom test.biom \\
            --metadata metadata.tsv \\
            --labels-col group \\
            --sample-id-col sample_id \\
            --embedding-birnn model.pt \\
            --plotfile-loss loss.png \\
            --plotfile-auc auc.png \\
            --glove-embedding result/embeddings_100.txt --d-model 100
    """
    set_log_level(kwargs['log'])
    params = {
        key: kwargs[key]
        for key in ('metadata', 'labels_col', 'sample_id_col', 'train_biom',
                    'valid_biom', 'test_biom', 'embedding_birnn',
                    'plotfile_loss', 'plotfile_auc', 'num_steps', 'p_drop',
                    'batch_size', 'd_model', 'd_ff', 'head_hidden', 'n_layers',
                    'n_heads', 'numb', 'lr', 'weight_decay', 'num_epochs',
                    'loss', 'alpha', 'glove_embedding', 'pred_out')
    }
    Attention_biom(**params)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
