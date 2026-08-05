# Trait analysis

This module contains trait association, prediction, and embedding/phylogeny
comparisons associated with Fig. 2A-D and Extended Data Fig. 5B,D-F. Figure
adoption and end-to-end run provenance remain unresolved.

| Subanalysis | Core files | Inputs | Outputs | Status |
| --- | --- | --- | --- | --- |
| Oxygen vector arithmetic | `oxygen_vector_arithmetic/Facultatively_Anaerobic_-_Aerobic_=_Anaerobic.ipynb` | trait labels, SNE, PhyloE, taxonomy under `data/` and sibling modules | retained cosine tables under `oxygen_vector_arithmetic/results/` | Legacy notebook; SNE/PhyloE roles require confirmation |
| PLS-DA | `plsda/run_plsda.R`, `plsda/plot_plsda.R` | BugBase, Traitar, BacDive labels and embeddings under `data/` | RDS model families remain under `data/`; figures belong under `plsda/results/` | Legacy source; required `traits_bugbase.csv` is absent |
| Trait prediction | `trait_prediction/run_vsearch.sh`, `Bugbase_predict.py`, `diff_datasets_predict.py` | FASTA, VSEARCH output, trait tables, taxonomy and embeddings | transferred traits and `auc_res.csv` under `trait_prediction/results/` | Legacy source; literal `xxx` paths remain |
| Phenotype prediction | `phenotype_prediction/phenotype.ipynb`, `plot_predict_res.R` | phenotype, trait, taxonomy and embedding tables | `predict_metabolics_res.csv` under `phenotype_prediction/results/` | Legacy source; notebook provenance unresolved |
| Embedding/phylogeny comparison | `embedding_phylogeny_comparison/plot_hexbin_embedding_cos_phylo.py` | shared SNE and retained Newick tree | figures belong under `embedding_phylogeny_comparison/results/` | Supporting legacy source |
| Cross-task plotting | `plot_traits.R` | tables and models from the analyses above | figures belong in the owning task's `results/figures/` | Ownership of individual panels is unresolved |

`data/` retains labels, model families, trait transfers and supporting trees.
The final SNE is consumed from `../sne_construction/data/`. Trait tables that
are produced by one task and plotted by another are not copied; consumers should
use explicit relative paths to the producer's result directory.

The module does not include a verified Traitar invocation, BacDive acquisition,
held-out trait-prediction folds, or adopted run records.
