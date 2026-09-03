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
        predict_bigg_gene.py : CarveMe DIAMOND step vs BiGG database
                                   │
                                   ▼
        per-genome BiGG gene tables  ──▶  observed trait tables
        (query_gene, BiGG_gene, score)     (mapping_bigg_gene_table.tsv,
                                             mapping_scores_table.tsv)
                                   │
                                   ▼
        picrust_predict.sh : PICRUSt2 hsp.py on reference tree bac.tre
                                   │
                                   ▼
        per-OTU predicted BiGG gene tables (data/OTU_bigg_gene/{otu}.tsv)
                                   │
                                   ▼
        build_metbolic_model_picrust2.py : CarveMe carve + gapfill on M3
        (gram-positive / gram-negative universe from data/fid_gram.json)
                                   │
                                   ▼
        per-OTU metabolic models (data/OTU_metabolic_model_M3/{otu}.xml)
                                   │
                                   ▼
        run_smetana.py : SMETANA global mode, M11 medium
                                   │
                                   ▼
        per-pair MIP / MRO tables (data/high_sim_res_M11.csv,
                                   data/low_SNE_res.csv)
                                   │
                                   ▼
        plot_resluts.R : group comparison and figures
                                   │
                                   ▼
        figures/metabolic_smetana.pdf
```

## Steps

### 1. BiGG prediction for reference genomes

`run_module_build.sh` annotates each reference genome with Prokka and builds
its reference CarveMe model (SLURM array job; genome list in
`genome_id.txt`). The resulting `.faa` proteomes feed the next script.

`predict_bigg_gene.py` runs CarveMe's DIAMOND step on the reference
proteomes and writes, per genome, a table of best bitscores over the full
BiGG gene universe (zero-filled when no hit):

```bash
python predict_bigg_gene.py genomes/ -o bigg_gene_scores -p 24   # protein FASTA dir
python predict_bigg_gene.py --self-test                          # unit test
```

Output row format (same as a PICRUSt2 trait table):

| query_gene | BiGG_gene | score |
| --- | --- | --- |
| gene_0 | STM_v1_0.STM0047 | 58 |
| gene_1 | STM_v1_0.STM0054 | 345 |

These per-genome tables are aggregated into the observed trait tables
`mapping_bigg_gene_table.tsv` (presence/absence) and
`mapping_scores_table.tsv` (scores) that the prediction step consumes.

### 2. BiGG prediction for OTU sequences

`picrust_predict.sh` uses the PICRUSt2 environment and its `hsp.py` hidden
state prediction algorithm to infer BiGG gene content and scores for OTUs
from the reference-genome trait tables and the reference tree
`data/bac.tre`:

```bash
source activate picrust2
hsp.py -t data/bac.tre --observed_trait_table mapping_bigg_gene_table.tsv \
    -o bigg_gene_predicted.tsv -p 4 -m mp -n
hsp.py -t data/bac.tre --observed_trait_table mapping_scores_table.tsv \
    -o scores_predicted.tsv -p 24 -m pic -n
```

The two predictions are combined into one table per OTU under
`data/OTU_bigg_gene/{otu}.tsv` with the same
`query_gene, BiGG_gene, score` columns as the reference tables.

### 3. Metabolic model construction for OTUs

`build_metbolic_model_picrust2.py` builds one CarveMe model per OTU from its
predicted BiGG gene table:

```bash
python build_metbolic_model_picrust2.py <otu_id>
```

For each OTU it reads `data/OTU_bigg_gene/{otu_id}.tsv`, scores the BiGG
reactions against CarveMe's GPR rules, picks the gram-positive or
gram-negative universe model according to `data/fid_gram.json`, carves the
model, and gap-fills it on the M3 medium defined in `data/media_db.tsv`.
The output SBML model is written to
`data/OTU_metabolic_model_M3/{otu_id}.xml`. The CPLEX solver is used.

`run_carveme.sh` dispatches this script over OTU id lists in
`data/split_files_add/` as a SLURM array job.

### 4. Pairwise metabolic interaction with SMETANA

`run_smetana.py` computes, for every pair in an input pair table, the
metabolic interaction potential (MIP) and metabolic resource overlap (MRO)
between the two OTU models:

```bash
python run_smetana.py <pair_table.csv> <run_dir>
```

For each row it copies the two models from `data/OTU_metabolic_model_M3/`
into a per-pair run directory, writes a SMETANA community file, and calls
SMETANA in `global` mode with the M11 medium from `data/media_db.tsv`,
excluding the inorganic compounds listed in `data/inorganic.txt`. Six pairs
are processed in parallel with joblib. Result files are written under
`data/smetana/results/` and reduced to the compact tables
`data/high_sim_res_M11.csv` (high-similarity pairs) and
`data/low_SNE_res.csv` (near-zero-similarity pairs), with columns
`mip, mro, cosine, phylo, co_occur`.

### 5. Visualization

`plot_resluts.R` produces `data/metabolic_smetana.pdf` from the two
compact result tables:

- panel a: boxplots comparing MRO and MIP between the cosine > 0.9 group
  and the cosine ≈ 0 group (t-test p-values annotated);
- panel b: within each group, scatter plots of MRO vs MIP split by the 100
  highest versus 100 lowest co-occurrence pairs, with marginal boxplots.

## Data layout

| Path | Content |
| --- | --- |
| `data/bac.tre` | reference phylogeny used for PICRUSt2-style prediction |
| `data/OTU_bigg_gene/{otu}.tsv` | predicted BiGG gene tables per OTU |
| `data/fid_gram.json` | gram status per OTU (chooses the CarveMe universe) |
| `data/media_db.tsv` | media library (M3 for gap-filling, M11 for SMETANA) |
| `data/inorganic.txt` | inorganic compounds excluded from SMETANA |
| `data/OTU_metabolic_model_M3/` | built OTU metabolic models (SBML) |
| `data/high_sim_res_M11.csv` | MIP/MRO results for high-similarity pairs |
| `data/low_SNE_res.csv` | MIP/MRO results for near-zero-similarity pairs |
| `figures/metabolic_smetana.pdf` | final figure |

## Dependencies

- Python: CarveMe, reframed, smetana, pandas, joblib, tqdm; DIAMOND on PATH;
  CPLEX solver
- conda environments: `picrust2` (step 2), `carvem_3.7` (steps 1 and 3)
- R: dplyr, ggplot2, ggpubr, aplot, cowplot, and the other packages loaded
  in `plot_resluts.R`

## Upstream / downstream context

Pair definitions and SNE similarities come from
[SNE construction](../sne_construction/README.md): the `cosine` and
`co_occur` columns in the result tables are SNE embedding cosine similarity
and ecological co-occurrence, respectively. The retained pair tables are the
adopted source data for the figure; their deterministic regeneration from
the upstream SNE and co-occurrence artifacts is not re-evidenced here.
