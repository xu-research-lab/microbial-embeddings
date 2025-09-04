# SNEs: Microbial Social Niches Learned from >210,000 Human Gut Microbiomes for Improve Deep Learning-based Disease Classification

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)`membed` package adapts Natural Language Processing techniques to create Social Niche Embeddings(SNEs) for microbes based on their co-occurrence patterns across samples. These embeddings provide ecological representations of microbial taxa based on their community context.

![SNE](img/img1.png)

## Human gut microbiome resource

* **Pre-training Microbiome Biom Table:** [pretraining_table_filter.biom](./data/pretraining_table_filter.biom)
  * **Description:** This BIOM-format file contains the raw input data for training. It is a large-scale feature table of Operational Taxonomic Unit (OTU) abundance, comprising **21,091 samples** and **14,093 distinct OTUs**. The data was derived from the processing of 16S rRNA gene sequences from a collection of diverse human gut microbiome samples. 
* **SNEs (100-dimensional):** [social_niche_embedding_100.txt](./data/social_niche_embedding_100.txt)
  * **Description:** This text file contains the pre-trained Social Niche Embeddings (SNEs). These embeddings were generated using the **GloVe** algorithm, which was trained on the microbial co-occurrence patterns within the BIOM table. Each row corresponds to a single OTU, followed by 100 floating-point values that constitute its learned vector. The file contains corresponding vectors for all 14,093 OTUs.


## Usage

**To reproduce our results, follow these steps to install and run the project.**

### Installation

We recommend using Conda to manage the environment and dependencies.

1. **Create or Update Conda Environment:** Use the provided file to create a new, clean environment:

   ```bash
   # Create a new environment named 'membed'
   conda env create --name membed --file requirements_dev.yml
   conda activate membed
   ```

   **Install the `membed` package:** Install in editable mode using pip (recommended for development):

   ```bash
   pip install -e .
   ```

   Alternatively, for a standard installation:

   ```bash
   pip install .
   ```

### Part 1: Generating SNEs

We adapt the **GloVe** (Global Vectors for Word Representation) model, a technique from natural language processing, to the field of microbiology. By learning from global co-occurrence statistics across more than 210,000 human gut microbiomes, the model embeds each microbe into a latent ecological space. 

- **Prerequisite:** A `biom` table (`table.biom`) containing OTU/ASV abundance data across many samples.

- **Step 1: Create a Feature Dictionary** This initial step creates a dictionary that maps each microbial feature (OTU/ASV) to a unique integer index and counts its non-zero occurrences across all samples. 

  ```bash
  membed dict -b table.biom -d feature-dict.csv
  ```

- **Step 2: Compute the Co-occurrence Matrix** A co-occurrence matrix is constructed by calculating the frequency or intensity of joint appearances for every pair of microbes across all samples.

  ```bash
  # Use the percentiled_co_abundance metric
  membed cooccur -b table.biom -c table.co --metric abundance_percentile --cpus 28
  ```

  Note: You can see all available metrics with `membed cooccur --help`. The `abundance_percentile`metric was selected as the primary method in this study. 

  + **Parameters**:
  + `-b`: **[Required]** The input path for the BIOM-format file.
    
  + `-c`: **[Required]** The output path for the co-occurrence matrix file.
    + `--metric`: Specifies the method used to quantify the association strength between microbes. Options include **binary** (presence/absence), **normalized abundance-based** (e.g., Bray-Curtis), and **percentile (rank)-based** metrics. The primary metric used in this study, `abundance_percentile`, is designed to overcome the limitations of other methods by assessing the **magnitude** of abundance for a more balanced evaluation, while also using a **percentile rank** transformation to ensure robustness against batch effects and extreme values.

- **Step 3: Calculate the `x_max` Hyperparameter **The `x_max` value is used to down-weight high-frequency co-occurrences during GloVe training, preventing them from dominating the loss function. For this purpose, the 80th percentile of the co-occurrence distribution is adopted as the threshold.

  ```bash
  membed build-x-max-file -c table.co -x xmax_file.npy --percentile_num 80
  ```

  + **Parameters**:
    + `-c`: **[Required]** The input co-occurrence matrix file generated in Step 2.
    + `-x`: **[Required]** The output path for the calculated `x_max` value, which will be saved in `.npy` format.
    + `--percentile_num`: The percentile of the co-occurrence distribution to use for calculating `x_max`. The recommended value is `80`.

- **Step 4: Train the GloVe Model to Generate Embeddings** The GloVe algorithm is trained using the co-occurrence matrix and the `x_max` value. 

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
  
  - **Parameters**:
    - `-d`: **[Required]** The input feature dictionary file from Step 1.
    - `-c`: **[Required]** The input co-occurrence matrix file from Step 2.
    - `-x`: **[Required]** The input `x_max` file from Step 3.
    - `-r`: **[Required]** The output directory where results will be stored. The SNEs file will be saved in this directory with the name `embeddings_{embedding-size}.txt`.
    - `--lr`: Learning Rate. This hyperparameter controls the step size of weight updates during model optimization.
    - `--embedding-size`: Embedding Dimensionality. This defines the length of the vector that represents each microbe. A dimensionality of `100` was found to perform well in our study.
    - `--iter`: Number of Iterations. 
    - `--cpus`: Number of CPU Cores. 

### Part 2: Downstream Classification using SNE

The **membed class-attention** module is an attention-based classification model based on SNEs.

- **Prerequisites:**

  - Training and testing data (`train.biom`, `test.biom`)
  - A metadata file (`metadata.tsv`) mapping sample IDs to their corresponding class labels.
  - Pre-trained SNEs: `embeddings_100.txt` 

- **Basic Example: Training an Attention-based Classifier (as used in our paper)**

  ```bash
  membed class-attention \
      -g embedding_100.txt \
      -tra_otu train.biom \
      -tes_otu test.biom \
      -m metadata.tsv \
      --labels_col group \
      --sample_id_col sample_id \
      -e attention_loo.pt \
      -ploss attention_loss.png \
      -pauc attention_auc.png \
      --num-steps 600 \
      --num-epochs 100 \
      --loss BCE_loss \
      --p-drop 0.4 \
      --d-ff 8 \
      --d-model 100 \
      --n-layers 1 \
      --n-heads 1 \
      --weight-decay 0.0001 \
      --lr 0.0001 \
      --batch-size 512 \
      --numb 8 \
  ```

+ **Argument Explanation**
  + `-g`, `-tra_otu`, `-tes_otu`, `-m`: Specify paths to the four mandatory input files: SNEs, training biom table, testing biom table, and metadata, respectively.
  + `-ploss`, `-pauc`, `-e`: Define the output paths for the loss curve plot, the AUC curve plot, and the model's attention weights.
  + `--labels_col group`: Informs the program that the column named `group` in the metadata file contains the classification labels.
  + `--num-epochs`, `--lr`, `--batch-size`, etc.: Set the model's hyperparameters, such as epochs, learning rate, batch size, and model dimensions, mirroring the definitions in your script.
  + `--d-model 100`: Defines the model's internal dimensionality. **This value must exactly match the dimension of the input embeddings.**

## Code Structure & Analysis Reproducibility

This repository is organized to reproduce every analysis presented in our paper. The `analysis/` directory contains subfolders, each corresponding to a specific figure or analytical theme.

- **`Pretraining_data_profile/`**: Scripts for building and profiling the pre-training dataset
- **`Co_occurence_method/`**: Scripts for the comparative analysis of co-occurrence metrics.
  - Comparative_analysis.R: A standalone R script for sampling OTUs from a single sample and generating network visualizations for each of the eight co-occurrence metrics.
- **`Simulation_experiments/`**: Scripts for validating the SNE framework using synthetic microbiome data.
- **`SNE_overview/`**: Code for visualizing the pre-trained Social Niche Embeddings
- **`Genome_collection_search/`**: Scripts for mapping OTUs to reference genomes
- **`Traits/`**: Scripts for analyzing the association between SNEs and microbial traits.
- **`Metabolic/`**: Scripts for metabolic interaction analysis using SMETANA.
- **`HGT/`**: Scripts for analyzing SNEs in relation to phylogeny, function, and Horizontal Gene Transfer.
- **`Disease_classification_loo/`**: Scripts for all disease and host phenotype classification experiments.
