"""Run CarveMe's DIAMOND step on genome proteomes and dump per-genome BiGG gene scores.

One CSV per genome, in the PICRUSt2 trait-table row format:

    query_gene,BiGG_gene,score
    gene_0,STM_v1_0.STM0047,58
    gene_1,STM_v1_0.STM0054,345
    gene_2,STM_v1_0.STM0059,0

Every BiGG gene in bigg_gprs.csv gets a row; no DIAMOND hit -> score 0.
No solver needed (nothing is optimized here), only `diamond` on PATH.

    python predict_bigg_gene.py genomes/ -o bigg_gene_scores -p 24
    python predict_bigg_gene.py --self-test
"""
import argparse
import os
from glob import glob

import pandas as pd
from carveme import config, project_dir
from carveme.reconstruction.diamond import run_blast, load_diamond_results

DIAMOND_ARGS = '--more-sensitive --top 10 --quiet'  # carveme defaults


def bigg_gene_universe():
    """All BiGG gene IDs, built the same way as carveme/reconstruction/scoring.py."""
    gprs = pd.read_csv(project_dir + config.get('generated', 'bigg_gprs'))
    return sorted({f"{m}.{g[2:]}" for m, g in zip(gprs.model, gprs.gene)})


def to_trait_rows(annotations, universe):
    """Best bitscore per BiGG gene, zero-filled over the full universe."""
    best = (annotations.sort_values('score', ascending=False)
                       .groupby('BiGG_gene')['score'].first())
    out = pd.DataFrame({'BiGG_gene': universe})
    out['score'] = out.BiGG_gene.map(best).fillna(0).astype(int)
    out.insert(0, 'query_gene', [f'gene_{i}' for i in range(len(out))])
    return out


def run_one(fasta, outdir, universe, dna, threads, diamond_args, force):
    name = os.path.splitext(os.path.basename(fasta))[0]
    blast_out = os.path.join(outdir, name + '.diamond.tsv')
    csv_out = os.path.join(outdir, name + '.csv')

    if os.path.exists(csv_out) and not force:
        return csv_out

    if not os.path.exists(blast_out) or force:
        code = run_blast(fasta, 'dna' if dna else 'protein', blast_out,
                         project_dir + config.get('generated', 'diamond_db'),
                         f'{diamond_args} -p {threads}', verbose=True)
        if code != 0:
            raise RuntimeError(f'diamond failed on {fasta} (exit code {code})')

    to_trait_rows(load_diamond_results(blast_out), universe).to_csv(csv_out, index=False)
    return csv_out


def self_test():
    ann = pd.DataFrame({'query_gene': ['q1', 'q2', 'q3'],
                        'BiGG_gene': ['M.a', 'M.a', 'M.b'],
                        'score': [100, 300, 50]})
    out = to_trait_rows(ann, ['M.a', 'M.b', 'M.c'])
    assert list(out.score) == [300, 50, 0], out           # best hit + zero fill
    assert list(out.query_gene) == ['gene_0', 'gene_1', 'gene_2'], out
    assert list(out.columns) == ['query_gene', 'BiGG_gene', 'score'], out
    print('self-test ok')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('input', nargs='?', help='protein FASTA, or a directory of them')
    p.add_argument('-o', '--outdir', default='bigg_gene_scores')
    p.add_argument('--dna', action='store_true', help='input is DNA FASTA (runs diamond blastx)')
    p.add_argument('--ext', default='.faa', help='file suffix in directory mode (default: .faa)')
    p.add_argument('-p', '--threads', type=int, default=4)
    p.add_argument('--diamond-args', default=DIAMOND_ARGS)
    p.add_argument('-f', '--force', action='store_true', help='recompute existing results')
    p.add_argument('--self-test', action='store_true')
    a = p.parse_args()

    if a.self_test:
        self_test()
        raise SystemExit

    if not a.input:
        p.error('need an input FASTA or directory')

    files = (sorted(glob(os.path.join(a.input, '*' + a.ext)))
             if os.path.isdir(a.input) else [a.input])
    if not files:
        p.error(f'no *{a.ext} under {a.input}')

    os.makedirs(a.outdir, exist_ok=True)
    universe = bigg_gene_universe()
    print(f'{len(universe)} BiGG genes x {len(files)} genomes -> {a.outdir}/')

    for i, f in enumerate(files, 1):
        print(f'[{i}/{len(files)}] {os.path.basename(f)}', flush=True)
        run_one(f, a.outdir, universe, a.dna, a.threads, a.diamond_args, a.force)

