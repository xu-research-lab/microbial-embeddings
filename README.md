# membed: Microbial Social Niche Embeddings from Large-Scale Microbiome Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)`membed` is a Python package that implements our pioneering **Social Niche Embedding (SNE)** framework. For a comprehensive description of the framework, please see our paper: [*“Microbial Social Niches Revealed from over 210,000 human gut microbiome samples”*].Grounded in the idea that a microbe's ecological role can be learned from its co-occurring neighbors, we adapted techniques from Natural Language Processing (NLP) to learn dense vector representations (i.e., "embeddings") for microbes. By training on a massive dataset of over 210,000 human gut microbiome samples, these embeddings position each microbe within a continuous "ecological space," revealing profound insights into their functions, interaction patterns, and connections to host health and disease.



## Key Findings 

Learned from an unprecedented amount of data, our SNEs are not just abstract vectors; they are rich with biological meaning:

- **🧬 Encodes Core Biological Traits**: SNEs effectively capture fundamental microbial phenotypes, including oxygen preference, Gram staining, cell shape, and key metabolic capabilities like sugar utilization.
- **🤝 Reveals Ecological Interactions**: SNE similarity is correlated with **competition** (Metabolic Resource Overlap) and **cooperation** (Metabolic Interaction Potential) between microbes. This allows for the prediction of potential microbial interactions directly from the embeddings.
- **🌳 Uncovers Insights Beyond Phylogeny**: SNEs capture ecological information that cannot be explained by evolutionary relatedness alone. SNEs are strongly associated with **Horizontal Gene Transfer (HGT)** events, especially for genes involved in environmental adaptation like antibiotic resistance and stress response.
- **🩺 Dramatically Improves Disease-State Classification**: Integrating pre-trained SNEs into machine learning models **substantially boosts diagnostic accuracy** for a wide range of diseases.
  - Outperforms traditional abundance-based models in diagnosing Inflammatory Bowel Disease (IBD) and Colorectal Cancer (CRC), achieving an overall AUC of 0.80 and 0.81, respectively.
  - Demonstrates remarkable generalizability, improving performance across 11 other microbiome-associated conditions, including Autism Spectrum Disorder, Parkinson's Disease, and Type 2 Diabetes.

## The `membed` Workflow

A typical workflow using the `membed` package involves two main stages:

1. **Generate SNEs**: Pre-train embeddings from your own large-scale microbiome dataset (e.g., a `biom` file with thousands of samples).
2. **Downstream Application**: Use the pre-trained embeddings as powerful features for various analyses, such as disease classification or biomarker discovery.

## Installation

We highly recommend using Conda to manage the environment and dependencies.

1. **Create or Update Conda Environment:** Use the provided file to create a new, clean environment:

   Bash

   ```
   # Create a new environment named 'membed'
   conda env create --name membed --file requirements_dev.yml
   conda activate membed
   ```

   Or, to update an existing environment:

   Bash

   ```
   conda activate your_existing_env_name
   conda env update -f requirements_dev.yml --prune
   ```

2. **Install the `membed` package:** Install in editable mode using pip (recommended for development):

   Bash

   ```
   pip install -e .
   ```

   Alternatively, for a standard installation:

   Bash

   ```
   pip install .
   ```



## Usage Tutorial

### Part 1: Generating Social Niche Embeddings (SNEs)

This workflow generates the core embedding vectors from a large-scale microbiome dataset.

- **Prerequisite:** A `biom` table (`table.biom`) containing OTU/ASV abundance data across many samples.

- **Step 1: Create a Feature Dictionary** This step creates a unique index for each OTU/ASV in your `biom` table, which is necessary for tracking features.

  ```bash
  membed dict -b table.biom -d feature-dict.csv
  ```

- **Step 2: Compute the Co-occurrence Matrix** This is a critical step where pairwise microbial co-occurrence is calculated. We recommend using the 

  `abundance_percentile` metric developed in our paper. 

  ```bash
  # Use the percentiled_co_abundance metric
  membed cooccur -b table.biom -c table.co --metric abundance_percentile
  ```

  Note: You can see all available metrics with `membed cooccur --help`. The 

  `abundance_percentile` metric is the core of our study. 

- **Step 3: Calculate the `x_max` Hyperparameter** The `x_max` value is used to down-weight high-frequency co-occurrences during GloVe training, preventing them from dominating the loss function. Our analysis identified the 80th percentile as a robust choice.

  ```bash
  membed build-x-max-file -c table.co -x xmax_file.npy --percentile_num 80
  ```

- **Step 4: Train the GloVe Model to Generate Embeddings** Use the co-occurrence matrix and `x_max` value to train the model and generate the final embeddings.

  Bash

  ```bash
  membed glove-train -d feature-dict.csv \
                     -c table.co \
                     -r ./result/ \
                     -x xmax_file.npy \
                     --lr 0.05 \
                     --embedding-size 100 \
                     --iter 100 \
                     --cpus 8
  ```

  - `-r ./result/`: Output directory where `embeddings_{embedding-size}.txt` will be saved.
  - `--embedding-size 100`: The dimensionality of the embedding vectors (100 was used in our paper).
  - `--iter 100`: The number of training iterations.

### Part 2: Using SNEs for Downstream Classification Tasks

Once you have the pre-trained `embeddings_100.txt`, you can leverage them as features in machine learning models.

- **Prerequisites:**

  - Training and testing data (`train.biom`, `test.biom`)
  - A metadata file (`metadata.tsv`) containing the sample labels (e.g., disease status)
  - Your pre-trained `embeddings_100.txt` file from Part 1

- **Example: Training an Attention-based Classifier (as used in our paper)**

  Bash

  ```bash
  membed class-attention \
      --train-biom train.biom \
      --test-biom test.biom \
      -m metadata.tsv \
      --labels_col group \
      --sample_id_col sample_id \
      -e embeddings_100.txt \
      -ploss attention_loss.png \
      -pauc attention_auc.png \
      --num-epochs 100 \
      --lr 0.0005 \
      --batch-size 64 \
      --embedding-size 100 \
      -g output_results/
  ```



## Code Structure & Analysis Reproducibility

This repository is organized to faithfully reproduce every analysis presented in our paper. The `analysis/` directory contains subfolders, each corresponding to a specific figure or analytical theme.

- **`Pretraining_data_profile/`**: Scripts for building and profiling the pre-training dataset, corresponding to **Figure S2**.
- **`Co_occurence_method/`**: Comparative analysis of the eight co-occurrence metrics, corresponding to **Figure S1**.
- **`Simulation_experiments/`**: Validation of the SNE framework using synthetic microbiome data, corresponding to **Figure 2**.
- **`SNE_overview/`**: Code for visualizing the pre-trained Social Niche Embeddings, corresponding to **Figure S4**.
- **`Genome_collection_search/`**: Scripts for mapping OTUs to reference genomes, supporting the analyses in **Figure S6**.
- **`Traits/`**: Analysis of the association between SNEs and microbial traits, corresponding to **Figure 3 & S5**.
- **`Metabolic/`**: Metabolic interaction analysis using SMETANA, corresponding to **Figure 4**.
- **`HGT/`**: Analysis of SNEs in relation to phylogeny, function, and Horizontal Gene Transfer, corresponding to **Figure 5**.
- **`Disease_classification_loo/`**: All disease and host phenotype classification experiments, corresponding to **Figure 6, 7, S8, & S9**.



## Data and Pre-trained Models

To facilitate reproducibility and community use, we provide:

- **Complete Source Code**: Available in this GitHub repository.
- **Pre-trained SNEs**: The 100-dimensional embeddings for all 14,067 OTUs are available.
- **Analysis Scripts**: All code required to reproduce our data pre-processing, embedding training, and downstream analyses.

All resources can be found in our GitHub repository: 

https://github.com/xu-research-lab/microbial-embeddings 



## Citation

If you use the `membed` package or our SNEs in your research, please cite our paper (placeholder, please update upon publication):

> Xu, Z. Z., et al. (2025). Microbial Social Niches Revealed from over 210,000 human gut microbiome samples. *Journal Name*, vol(issue), pages.

