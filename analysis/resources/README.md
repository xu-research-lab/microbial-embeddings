# Resources: pretraining data and genome-to-OTU mapping

This directory contains the pipelines that build the data resources of the
[SNEs](https://github.com/xu-research-lab/microbial-embeddings) project:

1. **16S database construction** ([`16S_database/`](16S_database/)) — turn raw
   SRA 16S amplicon reads into SILVA-mapped OTU tables, and assemble the
   pretraining BIOM from them.
2. **Reference-genome → OTU mapping** ([`genome_mapping/`](genome_mapping/)) —
   extract 16S sequences from high-quality reference genomes and map the SNE
   OTUs to genomes by sequence identity.

---

## 1. 16S database construction

The scripts process one study at a time; each study is a directory containing
its raw reads. The pipeline has three stages:

| Stage | Script | What it does |
|---|---|---|
| Download & convert | [`16S_database/process_data.sh`](16S_database/process_data.sh) | Downloads SRA accessions listed in `sample.txt` with `prefetch`, converts them to FASTQ with `fastq-dump --split-3` (20 jobs via GNU parallel), then launches the DADA2 step |
| DADA2 ASV inference | [`16S_database/process_forwards.R`](16S_database/process_forwards.R) | R script (`Rscript process_forwards.R <study_dir> trim`): quality filtering (`filterAndTrim`, `trimLeft = 24` for trimmed datasets), forward error model, `dada()`, ASV sequence table, chimera removal, and length filter (> 49 nt). Outputs `results/ASV.tsv`, `results/ASVs.fa`, `results/ASVs_counts.tsv`, `results/summary.tsv` |
| SILVA closed-reference clustering | [`16S_database/run_qiime_close_reference.sh`](16S_database/run_qiime_close_reference.sh) | QIIME 2 (`qiime2-amplicon-2024.2`): imports the ASV table and sequences, then `vsearch cluster-features-closed-reference` against SILVA 138.2 NR99 at 97% identity, and exports the clustered table |

### Pretraining dataset assembly

[`16S_database/get_pretraining_datasets.ipynb`](16S_database/get_pretraining_datasets.ipynb)
loads the merged compendium (`../../data/table_gut_all.biom`), removes the
samples used for the disease-classification analyses (matched against
`../../data/metadata_disease_classification.tsv`) so the pretraining set is
independent of the downstream classification evaluations, drops empty
features, and writes `../../data/gut_pretraining.biom`.

## 2. Reference-genome → OTU mapping

| Step | File | What it does |
|---|---|---|
| Genome collection | [`genome_mapping/Genome_collection.ipynb`](genome_mapping/Genome_collection.ipynb) | Filters the GTDB r220 metadata (`ar53_metadata_r220.tsv.gz`, `bac120_metadata_r220.tsv.gz`) to high-quality genomes: `checkm2_completeness > 90`, `checkm2_contamination < 5`, and not `derived from metagenome` → `data/high_quality_genome.tsv` |
| 16S extraction & vsearch | [`genome_mapping/run_vsearch.sh`](genome_mapping/run_vsearch.sh) | SLURM job: extracts 16S sequences from each genome with `barrnap`, concatenates them into `data/barrnap.fna`, then runs `vsearch --usearch_global` with the SNE OTU sequences (`data/feces_seq_16S_SLIVA.fasta`) as queries against the genome 16S database at 99% identity → `data/vsearch_blast.out` |
| Result QC | [`genome_mapping/genome_mapping_result.ipynb`](genome_mapping/genome_mapping_result.ipynb) | Parses the BLAST6 output, requires full query coverage (`q.start == 1`, coverage ≥ 1), keeps the best hit per genome, and produces the final mapping `data/vsearch_res.csv` |

The vsearch result tables are retained in two locations:
`genome_mapping/data/` (full inputs and outputs of this pipeline) and
`data/genome_mapping/` (the subset consumed by downstream analyses).

## 3. Data organization

```
16S_database/
├── process_data.sh                  # SRA download + fastq-dump + DADA2 launch
├── process_forwards.R               # DADA2: filtering, error model, ASV table
├── run_qiime_close_reference.sh     # QIIME 2 closed-reference clustering (SILVA 138.2)
└── get_pretraining_datasets.ipynb   # compendium -> gut_pretraining.biom (disease samples removed)

genome_mapping/
├── Genome_collection.ipynb          # GTDB r220 quality filtering -> high_quality_genome.tsv
├── run_vsearch.sh                   # barrnap 16S extraction + vsearch usearch_global (99% id)
├── genome_mapping_result.ipynb      # BLAST6 parsing, coverage QC -> vsearch_res.csv
└── data/                            # GTDB metadata, barrnap.fna, OTU query fasta,
                                     #   vsearch_blast.out, vsearch_res.csv

data/genome_mapping/                 # downstream-consumed subset of the mapping outputs
```

## 4. Environment

- R: `dada2`, `stringr` (R ≥ 4.x; the notebook steps use `biom`, `pandas`,
  `numpy` from the `membed` conda environment, see the [main
  README](../../README.md)).
- Shell: SRA Toolkit (`prefetch`, `fastq-dump`), GNU `parallel`.
- QIIME 2 environment `qiime2-amplicon-2024.2` with the `vsearch` plugin.
- `barrnap` and `vsearch` for genome 16S extraction and mapping.

## 5. How to reproduce

```bash
# Per-study 16S pipeline
bash 16S_database/process_data.sh          # edit study_id / file_path / path first
bash 16S_database/run_qiime_close_reference.sh

# Pretraining BIOM (after all studies are merged into table_gut_all.biom)
# run 16S_database/get_pretraining_datasets.ipynb  -> ../../data/gut_pretraining.biom

# Genome-to-OTU mapping
# run genome_mapping/Genome_collection.ipynb       -> data/high_quality_genome.tsv
sbatch genome_mapping/run_vsearch.sh
# run genome_mapping/genome_mapping_result.ipynb   -> data/vsearch_res.csv
```

The download, DADA2, and QIIME 2 scripts contain empty `study_id`, `file_path`,
and `path` variables and must be edited per study before running.
