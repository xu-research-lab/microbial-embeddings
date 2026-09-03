"""CPU smoke test for membed.otu_attention on the IBD toy cohort."""

import pytest

from membed.otu_attention import Attention_biom

pytestmark = pytest.mark.slow


def test_attention_cpu_smoke(data_dir, tmp_path):
    record, test_metrics = Attention_biom(
        metadata=str(data_dir / 'metadata_IBD.txt'),
        train_biom=str(data_dir / 'IBD_train.biom'),
        valid_biom=str(data_dir / 'IBD_valid.biom'),
        test_biom=str(data_dir / 'IBD_test.biom'),
        embedding_birnn=str(tmp_path / 'model.pt'),
        plotfile_loss=str(tmp_path / 'loss.png'),
        plotfile_auc=str(tmp_path / 'auc.png'),
        labels_col='group',
        sample_id_col='sample',
        num_steps=64,
        d_model=16,
        n_layers=1,
        batch_size=64,
        numb=-1,
        num_epochs=2,
        pred_out=str(tmp_path / 'pred'),
    )
    assert (tmp_path / 'model.pt').exists()
    assert (tmp_path / 'loss.png').exists()
    assert (tmp_path / 'auc.png').exists()
    assert (tmp_path / 'pred_valid.csv').exists()
    assert (tmp_path / 'pred_test.csv').exists()
    assert 0 <= test_metrics['auc'] <= 1
    assert 0 <= record['valid_auc'] <= 1
