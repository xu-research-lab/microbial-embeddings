# Metabolic interaction analysis

## Paper correspondence and purpose

Paper reference: Fig. 3A-B. This module contains legacy source for metabolic
resource overlap (MRO) and metabolic interaction potential (MIP) comparisons
for SNE-defined pairs using CarveMe/SMETANA candidates.

## Quick code map

| Scientific task | Location | Core/supporting code | Entry point | Status |
| --- | --- | --- | --- | --- |
| Metabolic model construction | `metabolic_model_construction/` | `build_metbolic_model_picrust2.py`, `run_carveme.sh` | UNKNOWN | legacy source; external inputs required |
| SMETANA interaction analysis | `smetana_interaction_analysis/` | `run_smetana.py`, `run_smetana.sh` | UNKNOWN | legacy source; external inputs required |
| Result visualization | `smetana_interaction_analysis/` | `plot_resluts.R` | UNKNOWN | legacy source; non-canonical result candidates |

## Subanalysis details

### Metabolic model construction

`build_metbolic_model_picrust2.py` reconstructs models with CarveMe, while
`run_carveme.sh` dispatches model identifiers. The source expects a CarveMe
installation and inputs under `data/`, including OTU-to-BiGG annotations,
Gram-status metadata, media definitions, and workload split files. Generated
model XML files are expected under the retained data layout. No executed-run
provenance or verified entry point is available.

### SMETANA interaction analysis and visualization

`run_smetana.py` prepares pair communities and calls SMETANA;
`run_smetana.sh` supplies pair-table and run-directory arguments.
`plot_resluts.R` performs the group comparisons and creates the MRO/MIP plots.
These sources require SMETANA, metabolic models, pair tables, media definitions,
and the inorganic-compound list retained or referenced under `data/`.

The plotting script and its two compact source-data candidates are kept in the
same subanalysis because visualization is the final stage of the SMETANA
analysis. Generated figures belong under
`smetana_interaction_analysis/results/figures/`.

## Results and provenance

The `smetana_interaction_analysis/results/` directory contains the source-data
candidates `high_sim_res_M11.csv` and `low_SNE_res.csv`. Their adoption,
deterministic generation, and executed-run provenance are not evidenced here.

Retained media, annotation, model, and raw SMETANA families are under `data/`.
Some pair tables and upstream identities referenced by the source are still
missing or external.

SNE context is documented in [SNE construction](../sne_construction/README.md).
No verified end-to-end entry point or rerun recipe is documented by this
module; task placement is scientific organization, not execution evidence.
