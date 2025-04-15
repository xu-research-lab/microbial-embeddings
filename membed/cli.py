# -*- coding: utf-8 -*-
"""Console script for membed."""
import sys
import os

from functools import partial

import functools
import wandb
import yaml

import click

from .glove import get_feature_dict, cooccur_workflow, train_glove_model
from .util import set_log_level
from .MLP import mlp
from .Attention_embedding import Attention_biom


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


jls_extract_var = 'Output the numbers of every feature.'


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
              help='Output co-occurrence matrix file path (format depends on metric).')
@click.option('-x',
              '--x-max-file',
              type=click.STRING,
              required=True,
              help='Output x_max value file path (used for GloVe training weight truncation).')
@click.option('--normalize', is_flag=True, default=False)
@click.option('--percentile', is_flag=True, default=False)
@click.option('--dense', is_flag=True, default=False,
              help='Use dense matrix computation (recommended for small datasets; sparse by default).')
@click.option('--metric', type=click.STRING, default='russell_rao',
                            help='''Similarity/distance metric (default: russell_rao):
1. **russell_rao** - Binary presence/absence similarity
   - Formula: (A∩B) / TotalSamples
   - Use Case: Basic co-occurrence frequency
   - Requires: Binary data (auto-converted)
   - Example: `--metric russell_rao`

2. **jaccard** - Jaccard similarity coefficient
   - Formula: (A∩B)/(A∪B)
   - Range: 0-1
   - Use Case: Presence/absence similarity
   - Requires: Binary data
   - Example: `--metric jaccard --dense`

3. **dice** - Dice/Sørensen-Dice coefficient
   - Formula: 2*(A∩B)/(|A|+|B|)
   - Range: 0-1
   - Use Case: Similarity weighting shared presences
   - Requires: Binary data
   - Example: `--metric dice --dense`

4. **faith** - Faith's phylogenetic similarity
   - Formula: (A∩B + (¬A∩¬B)/n)/n
   - Use Case: Phylogenetic community similarity
   - Requires: Binary data
   - Example: `--metric faith --dense`

5. **phi** - Phi correlation coefficient
   - Formula: (AD-BC)/√((A+B)(A+C)(B+D)(C+D))
   - Range: -1-1
   - Use Case: Binary correlation analysis
   - Requires: Binary data
   - Example: `--metric phi --dense`

6. **braycurtis** - Bray-Curtis dissimilarity
   - Formula: 1 - [2*∑min(u,v)]/(∑u+∑v)
   - Range: 0-1 (0=identical)
   - Use Case: Abundance-based dissimilarity
   - Requires: Normalized abundance data
   - Example: `--metric braycurtis --dense --normalize`

7. **abundance** - Abundance-weighted co-occurrence
   - Formula: ∑(min(u_i,v_i) - |u_i-v_i|*min(u_i,v_i))
   - Use Case: Weighted co-occurrence patterns
   - Requires: Raw abundance data
   - Example: `--metric abundance --dense`

8. **effect_size** - Co-occurrence effect size
   - Formula: (Observed - Expected)/MaxPossible
   - Use Case: Statistical over/under-cooccurrence
   - Requires: Hypergeometric expectation
   - Example: `--metric effect_size --dense`

9. **glove** - Position-weighted co-occurrence
   - Weight: 1/(position_distance)
   - Use Case: Embedding-focused relationships
   - Requires: Ordered abundance data
   - Example: `--metric glove --dense'''
)
@click.option('--cpus', type=click.INT, default=1)
@_common_parameters
def cooccur(**kwargs):
    """Compute pairwise cooccurrence of microbes in a biom table.

    Example:
    $ membed cooccur -b table.biom -c table.co -x xmax_file --metric russell_rao
    $ membed cooccur -b table.biom -c table.co -x xmax_file --metric russell_rao --dense
    $ membed cooccur -b table.biom -c table.co -x xmax_file --metric jaccard --dense
    $ membed cooccur -b table.biom -c table.co -x xmax_file --metric abundance --dense --percentile
    """
    set_log_level(kwargs['log'])
    params = {
        key: kwargs[key]
        for key in ('biom_file', 'cooccur_file', 'x_max_file', 'normalize',
                    'percentile', 'dense', 'metric', 'cpus')
    }
    cooccur_workflow(**params)


@main.command()
@click.option('-d',
              '--feature-dict',
              type=click.STRING,
              required=True,
              help='Intput the file of numbers of every feature.')
@click.option('-c',
              '--cooccur-file',
              type=click.STRING,
              required=True,
              help='file to store or load precomputed coccurrence.')
@click.option('-x',
              '--x-max-file',
              type=click.STRING,
              required=True,
              help='Output x-max for GloVe model.')
@click.option('-r',
              '--result',
              type=click.STRING,
              required=True,
              help='file to store embedding vector.')
@click.option('--lr', type=click.FLOAT, default=0.05)
@click.option('--embedding-size', type=click.INT, default=100)
@click.option('--iter', type=click.INT, default=50)
@click.option('--cpus', type=click.INT, default=4)
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
                     "x_max_file", "lr", "embedding_size", "iter",
                    "cpus")
    }
    train_glove_model(**params)


@main.command()
@click.option('-tra_otu',
              '--train-biom',
              type=click.STRING,
              required=True,
              help='file to store or load otu test.')
@click.option('-tes_otu',
              '--test-biom',
              type=click.STRING,
              required=True,
              help='file to store or load otu test.')
@click.option('-m',
              '--metadata',
              type=click.STRING,
              required=True,
              help='file to store or load glove embeding.')
@click.option('--labels_col', type=click.STRING, required=True)
@click.option('--sample_id_col', type=click.STRING, required=True)
@click.option('-e',
              '--embedding-birnn',
              type=click.STRING,
              required=True,
              help='file to store embedding_birnn.')
@click.option('-ploss',
              '--plotfile-loss',
              type=click.STRING,
              required=True,
              help='file to store loss plot.')
@click.option('-pauc',
              '--plotfile-auc',
              type=click.STRING,
              required=True,
              help='file to store auc plot.')
@click.option('--num-steps', type=click.INT, default=400)
@click.option('--p-drop', type=click.FLOAT, default=0.0)
@click.option('--d-ff', type=click.INT, default=64)
@click.option('--batch-size', type=click.INT, default=64)
@click.option('--d-model', type=click.INT, default=100)
@click.option('--n-layers', type=click.INT, default=2)
@click.option('--n-heads', type=click.INT, default=2)
@click.option('--numb', type=click.INT, default=1)
@click.option('--lr', type=click.FLOAT, default=0.0005)
@click.option('--weight-decay', type=click.FLOAT, default=0.0)
@click.option('--num-epochs', type=click.INT, default=1)
@click.option('--loss', type=click.STRING, default=None, help='balance loss.')
@click.option('--alpha', type=click.FLOAT, default=0.6, help='balance loss, alpha.')
@click.option('-g',
              '--glove-embedding',
              type=click.STRING,
              default=None,
              help='file to store or load glove embeding.')
@_common_parameters
def class_attention(**kwargs):
    """Run attention model.

    Example:
    """
    set_log_level(kwargs['log'])

    # with open(kwargs['config_file']) as f:
    #     sweep_config = yaml.load(f.read(), Loader=yaml.FullLoader)
    # sweep_id = wandb.sweep(sweep_config,
    #                        project=kwargs['wandb_project'],
    #                        entity="xu-lab")
    # wandb.agent(sweep_id,
    #             function=functools.partial(
    #                 Attention_biom, kwargs['config_file'], kwargs['metadata'],
    #                 kwargs['group'], kwargs['train_biom'], kwargs['test_biom'],
    #                 kwargs['embedding_birnn'], kwargs['num_steps'],
    #                 kwargs['p_drop'], kwargs['d_ff'],
    #                 kwargs['num_epochs'], kwargs['loss'],
    #                 kwargs['glove_embedding']))

    params = {
        key: kwargs[key]
        for key in ('metadata', 'labels_col','sample_id_col', 'train_biom', 'test_biom',
                    'embedding_birnn', 'plotfile_loss', 'plotfile_auc',
                    'num_steps', 'p_drop', 'd_ff', 'batch_size', 'd_model',
                    'n_layers', 'n_heads', 'numb', 'lr', 'weight_decay', 'num_epochs',
                    'loss', 'alpha', 'glove_embedding')
    }
    Attention_biom(**params)


@main.command()
@click.option('-tra_otu',
              '--train-biom',
              type=click.STRING,
              required=True,
              help='file to store or load otu test.')
@click.option('-tes_otu',
              '--test-biom',
              type=click.STRING,
              required=True,
              help='file to store or load otu test.')
@click.option('-m',
              '--metadata',
              type=click.STRING,
              required=True,
              help='file to store or load glove embeding.')
@click.option('--group', type=click.STRING, required=True)
@click.option('-ploss',
              '--plotfile-loss',
              type=click.STRING,
              required=True,
              help='file to store loss plot.')
@click.option('-pauc',
              '--plotfile-auc',
              type=click.STRING,
              required=True,
              help='file to store auc plot.')
@click.option('-mlpm',
              '--mlp-model',
              type=click.STRING,
              required=True,
              help='file to store model parameter.')
@click.option('--loss', type=click.STRING, default=None, help='balance loss.')
@click.option('-g',
              '--glove-embedding',
              type=click.STRING,
              default=None,
              help='file to store or load glove embeding.')
@click.option('--p-drop', type=click.FLOAT, default=0.2)
@click.option('--num-epochs', type=click.INT, default=400)
@click.option('--numb', type=click.INT, default=0)
@click.option('--weight-decay', type=click.FLOAT, default=0.001)
@click.option('--batch-size', type=click.INT, default=256)
@click.option('--num-hiddens', type=click.INT, default=64)
@click.option('--lr', type=click.FLOAT, default=0.0005)
@_common_parameters
def class_mlp(**kwargs):
    """Run BIRNN model.

    Example:
    """
    set_log_level(kwargs['log'])

    # with open(kwargs['config_file']) as f:
    #     sweep_config = yaml.load(f.read(), Loader=yaml.FullLoader)
    # sweep_id = wandb.sweep(sweep_config,
    #                        project=kwargs['wandb_project'],
    #                        entity="xu-lab")
    # wandb.agent(sweep_id,
    #             function=functools.partial(
    #                 mlp, kwargs['metadata'], kwargs['group'],
    #                 kwargs['train_biom'], kwargs['test_biom'],
    #                 kwargs['config_file'], kwargs['glove_embedding'],
    #                 kwargs['loss'], kwargs['num_epochs']))
    params = {
        key: kwargs[key]
        for key in ('metadata', 'group', 'train_biom', 'test_biom', 'numb',
                    'glove_embedding', 'plotfile_loss', 'plotfile_auc', 'mlp_model', 'p_drop',
                    'num_hiddens', 'batch_size', 'lr', 'weight_decay', 'num_epochs',
                    'loss')
    }
    mlp(**params)


if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
