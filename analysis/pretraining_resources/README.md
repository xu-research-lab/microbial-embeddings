# Pretraining resources

This module contains the retained material for construction and description of
the human gut microbiome compendium and OTU-to-genome mapping. It does not yet
provide a verified end-to-end preprocessing workflow.

| Subanalysis | Core files | Inputs | Outputs | Status |
| --- | --- | --- | --- | --- |
| Compendium overview | `compendium_overview/overview_data.ipynb`, `compendium_overview/gut_pretrain_overview.R` | `data/compendium/`, `data/profile/`, taxonomy mapping | profile tables under `data/profile/`; figures under `compendium_overview/results/figures/profile/` | Legacy source; adopted run unknown |
| Forward-read preprocessing | `profile_preprocessing/process_forwards.R` | A study working directory, `rum_sample.txt`, and raw reads | DADA2 intermediates and per-study `results/` in the supplied working directory | Legacy source; not a module-level rebuild entry point |
| Genome mapping | `genome_mapping/run_vsearch.sh` | `data/genome_mapping/` query FASTA and Barrnap database | VSEARCH output and mapping tables in `data/genome_mapping/` | Legacy source; representative-genome selection is not documented |

`data/compendium/` holds the retained BIOM and metadata inputs. `data/profile/`
holds overview tables and taxonomy data. `data/genome_mapping/` holds retained
mapping inputs and candidates. The final SNE and downstream modules consume
some of these retained files through explicit relative paths; they are not
duplicated by this README.

The included code does not establish the study manifest, cross-study QC,
closed-reference clustering, final filtering, or an adopted run record. Those
missing steps must be supplied before this module can reproduce the manuscript
compendium.
