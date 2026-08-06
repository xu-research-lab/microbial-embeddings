# Function, phylogeny, and HGT analysis

This module groups the retained function prediction, phylogeny comparison, HGT,
and Fig. 4 visualization material. No end-to-end entry point or adopted result
provenance is established.

| Subanalysis | Core files | Inputs | Outputs | Status |
| --- | --- | --- | --- | --- |
| Function prediction | `function_prediction/phylolm_subsystem_predict.R` | function table, phylogenetic tree and shared SNE | per-function RDS models under `data/phylolm/` | Legacy source; aggregation to `KO_phyloglm.csv` is absent |
| Phylogeny comparison | `phylogeny_comparison/run_SNE_phylo_function.py` | SNE, function profiles and tree | figures belong under `phylogeny_comparison/results/figures/` | Legacy source; dynamic function inputs are absent |
| HGT analysis | `hgt_analysis/HGT.ipynb`, `run_blastn.sh`, `run_HGT_predict.py` | genome pairs, BLAST resources, SNE, PhyloE and HGT labels | prediction candidate `data/hgt_predict_res_all.csv` | Legacy source; sample space and HGT aggregation are unresolved |
| Result visualization | `result_visualization/hgt.R`, `plot_results.R` | retained module data and compact HGT/taxonomy summaries | `result_visualization/results/tax_group.csv`, `result_visualization/results/hgt_plot_res.csv` | Legacy source; panel adoption unknown |

The HGT module uses the root-level canonical SNE input at
`../../data/social_niche_embedding_100.txt` (SHA-256
`374511e8610fa406d9d3e98f6dee74ce3485e379be64b7cf1e9516b02bd3a22f`),
resolved from the module root. This is the selected module input and does not
use the separate SNE-module copy. Function-model RDS files and intermediate
HGT prediction data remain in `data/`; compact tables used only for
visualization remain with `result_visualization/results/`.

The retained scripts are scientifically inconsistent: visualization code uses
Pearson-style statistics while earlier documentation describes Spearman-style
reporting; thresholds differ between scripts; and `plot_results.R` expects an
`SNE+PhyloE` group not emitted by `run_HGT_predict.py`. These are recorded
as unresolved evidence, not reconciled by reorganization.
