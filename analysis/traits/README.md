# Trait analysis

This directory contains all phenotype-related analyses.

## Directory structure

- `data/` — Data used for phenotype analysis and intermediate/generated data produced during the analysis.
- `traits_annotation/` — Scripts for phenotypic annotation of microbial OTUs.
  - `Bugbase_database.txt`, `Bugbase_predict.py`, `run_vsearch.sh`, `traitar_predict.sh`
  - Output files: `traits_predict_bacDive.csv`, `traits_predict_Bugbase.txt`, `traits_predict_Traitar.csv`
- `Facultatively_Anaerobic_-_Aerobic_=_Anaerobic.ipynb` — Natural-language-style analogy analysis (e.g. "Facultatively Anaerobic − Aerobic ≈ Anaerobic").
- `run_plsda.R` — PLS-DA analysis script.
- `traits_predict.ipynb` — Phenotype prediction notebook.
- `traits_results.ipynb` — Collection and visualization of all phenotype-related results.
