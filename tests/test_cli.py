"""Smoke tests for the membed command-line interface."""

import pytest
from click.testing import CliRunner

from membed.cli import main


@pytest.mark.parametrize('cmd', ['dict', 'cooccur', 'build-x-max-file',
                                 'glove-train', 'class-attention'])
def test_main_help(cmd):
    result = CliRunner().invoke(main, ['--help'])
    assert result.exit_code == 0
    assert 'Usage:' in result.output


@pytest.mark.parametrize('cmd', ['dict', 'cooccur', 'build-x-max-file',
                                 'glove-train', 'class-attention'])
def test_subcommand_help(cmd):
    result = CliRunner().invoke(main, [cmd, '--help'])
    assert result.exit_code == 0
    assert 'Usage:' in result.output


def test_cli_dict(data_dir, tmp_path):
    out = tmp_path / 'feature-dict.csv'
    result = CliRunner().invoke(
        main, ['dict', '-b', str(data_dir / 'test_raw.biom'), '-d', str(out)])
    assert result.exit_code == 0, result.output
    assert out.read_text().splitlines() == [
        'f0 2', 'f1 2', 'f2 3', 'f3 2', 'f4 3']


def test_cli_cooccur(data_dir, tmp_path):
    out = tmp_path / 'table.co'
    result = CliRunner().invoke(
        main, ['cooccur', '-b', str(data_dir / 'test_raw.biom'),
               '-c', str(out), '--metric', 'russell_rao', '--cpus', '1'])
    assert result.exit_code == 0, result.output
    assert out.stat().st_size > 0
