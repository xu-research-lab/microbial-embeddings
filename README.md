# membed: Global Microbial Embedding from Large-Scale Microbiome Data Reveals the Gut Microbiome Assembly Rules beyond Phylogeny

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) `membed` is a Python package implementing the **GloME (Global Microbial Embedding)** framework described in the paper "[*GloME: Global Microbial Embedding on Large-Scale Microbiome Datasets Unmasks the Gut Microbiome Assembly Rules beyond Phylogeny*]". It leverages large-scale microbiome datasets to generate meaningful vector representations (embeddings) of microbial taxa (OTUs/ASVs) based on their ecological co-occurrence patterns.

These embeddings capture complex ecological relationships and have been shown to enhance downstream analyses such as disease classification, biomarker identification, and understanding microbial community assembly rules, including the influence of Horizontal Gene Transfer (HGT).

## Key Features

* **GloME Embedding Generation:** Implements the core GloME algorithm (adapted GloVe) using a novel abundance-percentile co-occurrence metric suitable for microbiome data.
* **Large-Scale Data Handling:** Designed to process extensive microbiome datasets (e.g., `biom` format tables).
* **Downstream Model Integration:** Includes tools and examples for using embeddings in machine learning models (Attention-based Deep Learning, MLP).
* **Ecological Insights:** Facilitates analyses related to microbial co-occurrence networks, niche overlap, and distinguishing health/disease states.

## Installation

It is recommended to manage dependencies using Conda.

1.  **Create/Update Conda Environment:**
    You can create a new environment using the provided development requirements file:
    ```bash
    # Create a new environment named 'membed' (or choose your own name)
    conda env create --name membed --file requirements_dev.yml
    conda activate membed
    ```
    Or, if you have an existing environment, update it:
    ```bash
    conda activate your_existing_env_name
    conda env update -f requirements_dev.yml --prune
    ```

2.  **Install `membed`:**
    Install the package using pip (editable mode recommended for development):
    ```bash
    pip install -e .
    ```
    Alternatively, for a standard installation:
    ```bash
    python setup.py install
    # or
    # pip install .
    ```

## Workflow & Usage

### 1. Generating GloME Embeddings

This workflow generates the core co-occurrence-based embeddings from a large microbiome dataset (represented as a `biom` table).

* **Prerequisite:** Your input `biom` table (`table.biom`) should contain OTU/ASV abundance data across many samples. Feature IDs should ideally be consistent (e.g., mapped to a reference database like Greengenes or SILVA if comparing across datasets later).

* **Step 1: Create Feature Dictionary:** (Optional but recommended for tracking) Generate a mapping of feature IDs to indices.
    
    ```bash
    membed dict -b table.biom -d feature-dict.csv
    ```
    
* **Step 2: Calculate Co-occurrence Matrix:** Compute the pairwise microbial co-occurrence matrix using the specified metric (e.g., `russell_rao`, or the manuscript's `abundance_percentile` if implemented). This step also determines the `xmax` scaling factor for GloVe.
    
    ```bash
    # Example using russell_rao metric (adjust --metric as needed)
    membed cooccur -b table.biom -c table.co -x xmax_file.npy --metric russell_rao
    ```
    *Note: Refer to `membed cooccur --help` for available metrics and options. The manuscript primarily uses the "abundance-percentile" method.*
    
* **Step 3: Train GloVe Model:** Train the GloVe model using the co-occurrence matrix to generate the final embedding vectors.
    
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
    * `-r ./result/`: Output directory for embeddings (`embeddings.txt` or similar).
    * `--embedding-size`: Dimensionality of the embedding vectors (e.g., 100 as used in the paper).
    * `--iter`: Number of training iterations.
    * Adjust parameters (`--lr`, `--embedding-size`, `--iter`, `--cpus`) as needed.

### 2. Using Embeddings for Downstream Classification

Once you have the pre-trained embeddings (`embeddings.txt`), you can use them as input features for downstream machine learning tasks, such as disease classification.

* **Prerequisites:**
    
    * Training data (`train.biom`, `train_metadata.csv`)
    * Testing data (`test.biom`, `test_metadata.csv`)
    * Pre-trained embedding file (`embeddings.txt` from Step 1).
    * Metadata file should contain the target variable (e.g., disease status) specified by the `--group` parameter.
    
* **Example: Training an Attention-based Classifier:**
    
    ```bash
    membed class-attention \
        --train-biom path/to/train.biom \
        --test-biom path/to/test.biom \
        -m path/to/metadata.csv \
        --group <column_name_for_classification_target> \
        -e path/to/embeddings.txt \
        -ploss path/to/attention_loss.png \
        -pauc path/to/attention_auc.png \
        -g path/to/output_results/ # Output directory
    ```
    
* **Example: Training an MLP Classifier:**
    
    ```bash
    membed class-mlp \
        --train-biom path/to/train.biom \
        --test-biom path/to/test.biom \
        -m path/to/metadata.csv \
        --group <column_name_for_classification_target> \
        -e path/to/embeddings.txt \
        # Add other MLP-specific parameters (e.g., output paths) - check help
        # membed class-mlp --help
    ```

## Applications

The generated embeddings can be used for various microbiome analyses:

* **Disease Classification:** Improve prediction accuracy for conditions like IBD, CRC, etc.
* **Biomarker Discovery:** Identify microbes or microbial ecological roles associated with specific phenotypes.
* **Microbial Co-occurrence Network Analysis:** Explore ecological interactions.
* **Niche Analysis:** Understand niche overlap and competition/cooperation patterns.
* **HGT Prediction:** Investigate potential horizontal gene transfer events correlated with ecological similarity.
* **Community Assembly:** Study the rules governing how microbial communities form.

## For Developers

* **Environment Setup:** Use the `requirements_dev.yml` file with Conda as described in the Installation section.
* **Running Tests (pytest):**
    
    ```bash
    # Run all tests with verbose output
    pytest -vv -rA --doctest-modules --doctest-continue-on-failure
    
    # Run a specific test function
    pytest -vv --log-cli-level debug -k <test_function_name>
    
    # Enable tqdm progress bars during tests (useful for long tests)
    pytest -s -vv ...
    ```
* **Profiling (CPU and Memory):**
    Example command for profiling (adjust paths and parameters):
    
    ```bash
    pytest --log-level debug -s -k <profiling_test_function> \
           --n_jobs 3 \
           --table tests/data/example_table.biom \
           --metric binary \
           --mpo /tmp/mem.log \
           --lpo /tmp/time.log
    ```

## TODO / Future Work

* [x] 

If you use `membed` or the GloME framework in your research, please cite the original paper: