# Metabolic interaction analysis

## What we do here

This module computes pairwise metabolic interactions between microbial OTUs.
Because OTUs lack genome-scale metabolic models, we first **predict** a BiGG
gene profile for every OTU (from reference genomes, via a PICRUSt2-style
ancestral-state reconstruction), then **build** a CarveMe metabolic model per
OTU from that profile, and finally run **SMETANA** on pairs of models to score
their metabolic interaction potential (MIP, cooperation) and metabolic
resource overlap (MRO, competition).

The scores are then compared between SNE-defined pair groups (paper Fig. 3A-B):
pairs with high embedding similarity (cosine > 0.9) versus pairs with
near-zero similarity (cosine ≈ 0), and within each group, high versus low
co-occurrence pairs.

## Pipeline overview

```
reference genomes ──prokka──▶ proteomes (.faa)
                                   │
                                   ▼
        DIAMOND vs BiGG ─▶ per-genome hit tables (data/blast_output_bigg/)
                                   │
                                   ▼
        build_bigg_gene_tables.ipynb, part 1
          ─▶ genome x gene trait tables (data/mapping_bigg_gene_table.tsv,
                                          data/mapping_scores_table.tsv)
                                   │
                                   ▼
        picrust_predict.sh : PICRUSt2 hsp.py on reference tree data/bac.tre
          ─▶ data/bigg_gene_predicted.tsv, data/scores_predicted.tsv
                                   │
                                   ▼
        build_bigg_gene_tables.ipynb, part 2
          ─▶ per-OTU BiGG gene tables (data/OTU_bigg_gene/{otu}.tsv)
                                   │
                                   ▼
        build_metbolic_model_picrust2.py : CarveMe carve + gapfill on M3
          ─▶ per-OTU metabolic models (data/OTU_metabolic_model_M3/{otu}.xml)
                                   │
                                   ▼
        run_smetana.py : SMETANA global mode, M11 medium
          ─▶ per-pair results, reduced to data/high_sim_res_M11.csv and
             data/low_SNE_res.csv
                                   │
                                   ▼
        plot_resluts.R ─▶ data/metabolic_smetana.pdf
```

Every intermediate is provided (in the repository or in the Release), so any
step can be run on its own. Only step 5 is needed to redraw the figure.

## Steps

Run everything from this directory.

### 0. Fetch the large inputs

Six inputs exceed GitHub's file size limit and are published as the
[`metabolic-interaction-data-v1`](https://github.com/xu-research-lab/microbial-embeddings/releases/tag/metabolic-interaction-data-v1)
GitHub Release instead:

| Under `data/` | Size | Needed for |
| --- | --- | --- |
| `blast_output_bigg/` | 853 MB packed | step 1 |
| `mapping_bigg_gene_table.tsv`, `mapping_scores_table.tsv` | 2.8 GB, 3.0 GB | step 2 |
| `bigg_gene_predicted.tsv`, `scores_predicted.tsv` | 728 MB, 793 MB | step 2 |
| `OTU_metabolic_model_M3/` | 3.4 GB packed | step 4 |

```bash
bash ../../download_release_data.sh metabolic_interaction
```

The script (at the repository root) downloads, checksums and unpacks everything into `data/`. The model
archive is split into two parts to fit the Release size limit; the script joins
them. Finished files are skipped, so it is safe to re-run. Step 5 (the figure)
needs none of these files.

### 1. BiGG gene tables for the reference genomes

The reference genomes were annotated with Prokka beforehand, and their
proteomes were searched against the BiGG gene database with DIAMOND. Neither
the proteomes nor those job scripts are part of this repository. The results
are the per-genome hit tables `data/blast_output_bigg/{genome}.tsv`
(columns `query_gene, BiGG_gene, score`; one row per BiGG gene with a hit,
`score` is the DIAMOND bitscore).

The first part of `build_bigg_gene_tables.ipynb` merges these tables into the
two genome x BiGG gene trait tables that PICRUSt2 reads:

- `data/mapping_bigg_gene_table.tsv`: 1 if the genome has a hit for the gene, else 0;
- `data/mapping_scores_table.tsv`: the bitscore, else 0.

`data/genome_id.txt` lists 26,868 genomes; 13 of them have no hit table and
are skipped (the notebook prints how many), so both tables have 26,855 rows.

`predict_bigg_gene.py` redoes the DIAMOND step with CarveMe's own settings:

```bash
python predict_bigg_gene.py genomes/ -o bigg_gene_scores -p 24   # protein FASTA dir
python predict_bigg_gene.py --self-test                          # unit test
```

Its output is **not** a drop-in replacement for `data/blast_output_bigg/`. It
writes one `.csv` per genome with a row for every BiGG gene (score 0 when there
is no hit) and renames the query genes `gene_0, gene_1, ...`. The notebook
counts every row as a hit, so this output would mark every gene as present;
drop the zero-score rows first.

### 2. BiGG gene tables for the OTUs

`picrust_predict.sh` runs PICRUSt2's `hsp.py` on the reference tree
`data/bac.tre` to predict, for each OTU placed on the tree, gene presence
(maximum parsimony) and bitscore (phylogenetic independent contrasts):

```bash
bash picrust_predict.sh   # needs the picrust2 conda environment
```

It writes `data/bigg_gene_predicted.tsv` and `data/scores_predicted.tsv` (one
row per OTU, plus the `metadata_NSTI` and `closest_reference_genome` columns).

The second part of `build_bigg_gene_tables.ipynb` turns them into one table per
OTU, `data/OTU_bigg_gene/{otu}.tsv`. It reads the SNE vocabulary from
`../../data/social_niche_embedding_100.txt` and:

1. restores the OTU IDs. In `data/bac.tre`, and therefore in the `hsp.py`
   outputs, every `U` of an OTU ID was written as `T` (`ABOU02000049...`
   became `ABOT02000049...`). This affects 2,498 of the OTUs kept below; the
   notebook maps the IDs back to the SNE IDs;
2. keeps the OTUs that are in the SNE vocabulary and have NSTI < 2 (14,039
   OTUs);
3. for each OTU keeps the genes that are predicted present **and** have a
   predicted score above 0.

This reproduces the committed tables byte for byte. Their format is the one
CarveMe reads, with two quirks that do not affect CarveMe: the files are
comma-separated despite the `.tsv` suffix and start with an unnamed row-index
column, and the `gene_N` numbers have gaps where genes were dropped.

### 3. Metabolic models for the OTUs

`build_metbolic_model_picrust2.py` builds one CarveMe model per OTU from its
table:

```bash
python build_metbolic_model_picrust2.py <otu_id>
```

It scores the BiGG reactions against CarveMe's GPR rules, carves the model and
gap-fills it on the M3 medium from `data/media_db.tsv`, using the CPLEX solver.
The universe model comes from `data/fid_gram.json`: gram-positive (3,068 OTUs)
or gram-negative (887 OTUs). The other 10,084 OTUs are not listed there and use
CarveMe's default bacterial universe. The model is written to
`data/OTU_metabolic_model_M3/{otu_id}.xml`.

The script builds one OTU per call; the 14,039 OTUs were run as a SLURM array
job whose submission script is not part of this repository. 13,942 models were
built; the other 97 OTUs have none. The models are in the Release (step 0).

### 4. Pairwise interactions with SMETANA

```bash
python run_smetana.py <pair_table.csv> <run_dir>
```

The pair table is a CSV with a header row and exactly two columns, one OTU ID
each. For every pair the script copies the two models into
`<run_dir>/<row>/`, writes a SMETANA community file and runs SMETANA in
`global` mode on the M11 medium from `data/media_db.tsv`, ignoring the
compounds in `data/inorganic.txt`. Six pairs run in parallel. Pairs involving
one of the 97 OTUs without a model are skipped with a message. Results go to
`data/smetana/results/{otu1}_{otu2}_M11_output_global.tsv`.

The pair tables, and the script that reduced the SMETANA results to the two
tables used for the figure, are not in this repository. The two tables are
provided as-is:

| File | Pairs | Columns |
| --- | --- | --- |
| `data/high_sim_res_M11.csv` | 67,425 pairs with SNE cosine from 0.60 to 0.995; the figure uses the 1,099 with cosine > 0.9 | `mip, mro, cosine, phylo, co_occur` |
| `data/low_SNE_res.csv` | 17,399 pairs with SNE cosine between -0.001 and 0.001 | `mip, mro, co_occur, cosine` |

### 5. Figure

```bash
LC_ALL=en_US.UTF-8 Rscript plot_resluts.R
```

The script writes `data/metabolic_smetana.pdf`. The group labels contain "≈", so it
needs a UTF-8 locale (it stops with a message otherwise) and an R build with
cairo support (`capabilities("cairo")`).

- Panel a: MRO and MIP of the cosine > 0.9 group versus the cosine ≈ 0 group,
  with Welch t-test p-values. The y-axes are zoomed; no points are dropped
  from the boxplots.
- Panel b: within each group, MRO versus MIP for the 100 pairs with the highest
  and the 100 with the lowest co-occurrence, with marginal boxplots. MIP values
  are jittered by up to ±0.2 for display (fixed seed).

In the cosine ≈ 0 group, 12,693 of the 17,399 pairs have `co_occur` = 0, so
its "100 lowest" pairs are simply the first 100 of those rows in file order.

## Data layout

| Path | Content |
| --- | --- |
| `data/genome_id.txt` | reference genomes read by the notebook |
| `data/blast_output_bigg/{genome}.tsv` | DIAMOND hits against BiGG per reference genome (Release) |
| `data/mapping_bigg_gene_table.tsv`, `data/mapping_scores_table.tsv` | genome x BiGG gene presence and bitscore tables (Release) |
| `data/bac.tre` | tree of the 26,868 reference genomes with the OTUs placed on it (OTU IDs with `U` written as `T`), used by `hsp.py` |
| `data/bigg_gene_predicted.tsv`, `data/scores_predicted.tsv` | `hsp.py` predictions per OTU (Release) |
| `data/OTU_bigg_gene/{otu}.tsv` | BiGG gene table per OTU, 14,039 files |
| `data/fid_gram.json` | gram status for 3,955 OTUs (chooses the CarveMe universe) |
| `data/media_db.tsv` | media library (M3 for gap-filling, M11 for SMETANA) |
| `data/inorganic.txt` | compounds ignored by SMETANA (`o2`, `h`, `h2o`, `pi`) |
| `data/OTU_metabolic_model_M3/{otu}.xml` | OTU metabolic models (SBML), 13,942 files (Release) |
| `data/high_sim_res_M11.csv`, `data/low_SNE_res.csv` | MIP/MRO results used for the figure |
| `data/metabolic_smetana.pdf` | the figure, written by step 5 |

## Dependencies

- Step 1 notebook part and step 2 notebook part: `pandas`, `numpy`, `tqdm`.
- Step 1 (`predict_bigg_gene.py`) and step 3: CarveMe with DIAMOND on `PATH`,
  `reframed`, the CPLEX solver (conda environment `carvem_3.7`).
- Step 2 (`picrust_predict.sh`): conda environment `picrust2`.
- Step 4: `smetana`, `reframed`, CPLEX, `pandas`, `joblib`, `tqdm`.
- Step 5: R with cairo support and the packages `dplyr`, `ggplot2`, `cowplot`,
  `aplot`, `ggplotify`, `RColorBrewer`.

## Upstream / downstream context

Pair definitions and SNE similarities come from
[SNE construction](../sne_construction/README.md): the `cosine` and
`co_occur` columns in the result tables are SNE embedding cosine similarity
and ecological co-occurrence, respectively. The pair tables are provided
as-is: the scripts that generated them from the upstream SNE and
co-occurrence outputs are not part of this repository.
