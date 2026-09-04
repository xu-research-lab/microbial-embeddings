#!/usr/bin/env python
import numpy as np
import biom
from pathlib import Path
from joblib import Parallel, delayed

SIZES  = [1000, 2000, 5000, 10000, 20000, 40000, 80000, 160000]
LABELS = ["1k", "2k", "5k", "1w", "2w", "4w", "8w", "16w"]

def process_n(n, table_path, data_dir):
    rng = np.random.default_rng(seed=n)          # reproducible per replicate
    table = biom.load_table(table_path)          # load inside worker, avoid pickling
    sid = table.ids(axis="sample")

    out_dir = Path(data_dir) / f"data_{n}"
    out_dir.mkdir(parents=True, exist_ok=True)

    order = rng.permutation(sid)                 # shuffle once; prefixes are nested subsets

    for size, label in zip(SIZES, LABELS):
        if size > len(sid):
            print(f"[rep {n}] skip {label}: only {len(sid)} samples available")
            continue

        agsub = table.filter(order[:size], axis='sample', inplace=False)
        agsub.remove_empty(axis='observation', inplace=True)

        prevalence = agsub.nonzero_counts(axis='observation')
        fid = agsub.ids(axis='observation')
        keep = fid[prevalence >= 100]   # relative prevalence threshold
        agsub.filter(keep, axis='observation', inplace=True)
        agsub.remove_empty()

        outpath = out_dir / f"subset_table_{label}.biom"
        with biom.util.biom_open(str(outpath), 'w') as out:
            agsub.to_hdf5(out, f"subsampled from gut_pretraining (rep {n}, {label})")

        print(f"[rep {n}] {label}: {agsub.shape[1]} samples x {agsub.shape[0]} features")

if __name__ == '__main__':
    data_dir = "Data/pretraining_datasize/subset/"
    table_path = "../../data/gut_pretraining.biom"
    Parallel(n_jobs=5)(
        delayed(process_n)(n, table_path, data_dir) for n in range(1, 6)
    )