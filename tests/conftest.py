from pathlib import Path
import logging

import pytest

from membed.cooccur_embedding import read_biom

logger = logging.getLogger(__name__)


@pytest.fixture(scope='package')
def data_dir(request):
    '''Dir where test data are stored.'''
    logger.debug('setup fixture: data_dir')
    yield Path(request.fspath.dirname) / 'data'
    logger.debug('teardown fixture: data_dir')


@pytest.fixture
def biom_table(data_dir):
    def _biom_table(table_name):
        return read_biom(data_dir / table_name)
    yield _biom_table


def pytest_configure(config):
    config.addinivalue_line(
        'markers', 'slow: long-running tests (e.g. CPU training smoke tests)')
