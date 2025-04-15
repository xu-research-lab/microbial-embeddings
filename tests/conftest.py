from pathlib import Path
import logging

import pytest
from sklearn.metrics.pairwise import _VALID_METRICS

from membed.glove import read_biom

logger = logging.getLogger(__name__)

# persist and tear down only after running all tests
@pytest.fixture(scope='package')
def data_dir(request):
    '''Dir where test data are stored.'''
    logger.debug('setup fixture: data_dir')
    yield Path(request.fspath.dirname) / 'data'
    logger.debug('teardown fixture: data_dir')


@pytest.fixture
def biom_table(data_dir):
    def _biom_table(table_name, **kwargs):
        return read_biom(data_dir / table_name, **kwargs)
    yield _biom_table


def pytest_addoption(parser):
    parser.addoption(
        "--table",
        metavar='FILE',
        action="store",
        #required=True,
        help="input biom table",
        type=Path)
    parser.addoption(
        "--dense",
        action="store_true",
        help="convert biom table to dense array?")
    parser.addoption(
        "--cpus",
        metavar='N',
        action="store",
        default=1,
        help="number of CPU cores to use (default 1)",
        type=int)
    parser.addoption(
        "--metric",
        metavar='METRIC',
        action="store",
        default='abundance',
        choices=['abundance', 'binary'] + _VALID_METRICS,
        help="cooccur metric",
        type=str)
    parser.addoption(
        "--lpo",
        metavar='FILE',
        action="store",
        default=None,
        help="line_profiler output file path",
        type=Path)
    parser.addoption(
        "--mpo",
        metavar='FILE',
        action="store",
        default=None,
        help="memory_profiler output file path",
        type=Path)

@pytest.fixture
def table(request):
    return request.config.getoption("--table")

@pytest.fixture
def dense(request):
    return request.config.getoption("--dense")

@pytest.fixture
def metric(request):
    return request.config.getoption("--metric")

@pytest.fixture
def cpus(request):
    return request.config.getoption("--cpus")

@pytest.fixture
def mpo(request):
    return request.config.getoption("--mpo")

@pytest.fixture
def lpo(request):
    return request.config.getoption("--lpo")
