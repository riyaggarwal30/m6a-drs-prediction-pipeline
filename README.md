# Data Processing, Alignment, and Ground-Truth Curation

This section details the preprocessing pipeline used to combine sequencing runs, normalize coordinates, integrate orthogonal GLORI ground-truth datasets, and handle transcriptomic coverage biases to curate a high-confidence dataset for model training.

## 1. BAM Merging and Unified Pileup

To maximize read depth and preserve coverage across transcripts, raw BAM files from all four flow cells were concatenated using `samtools merge`, ensuring all overlapping reads across runs are evaluated together per site.

- **Script:** `scripts/01_merge_and_pileup.sh`
- **Output:** `data/processed/hsap_modkit_pileup.bed`

## 2. Header Normalization and Coverage Filtering

To ensure compatibility with standard reference annotations and GLORI coordinates, bedMethyl chromosome headers were normalized (converting `HSAP_` prefixes and aligning mitochondrial naming). Next, the file was filtered to retain high-confidence sites with a read coverage of \(\geq 10\times\) (column 10).

- **Script:** `scripts/02_filter_coverage.sh`
- **Output:** `data/processed/hsap_modkit_pileup_cov10.bed`

## 3.  Pipeline Execution via Jupyter Notebooks
All data processing, coordinate harmonization, feature engineering, and matrix exports are consolidated into two sequential Jupyter notebooks located in the scripts/ directory:


### 1. Training Set Pipeline (`scripts/01_build_train_set.ipynb`)

This notebook executes the complete end-to-end training data preparation workflow:

- **Coordinate Harmonization:** Standardizes DRS bedMethyl files (converting 0-based half-open coordinates to 1-based indexing) to match GLORI reference coordinates.
- **GLORI Union & Intersection:** Standardizes and merges positive and negative GLORI replicates, verifies zero conflicts, and performs an inner join with filtered DRS coverage data.
- **Dataset Balancing:** Downsamples the high-confidence intersection to create a balanced training pool (`balanced_ground_truth_4k.tsv`).
- **Feature Extraction & Aggregation:** Merges per-read modification probability statistics (`mod_prob_median`, `sd`, `q25`, `q75`) and Phred base quality scores (`base_qual_q25`, `q50`, `q75`) grouped by genomic locus.
- **Multi-Model Export:** Generates and exports the clean training matrices for all three feature ablation configurations (Model 1: QC/Sequence only; Model 2: Probabilities baseline; Model 3: Comprehensive feature set excluding modified fraction).

### 2. Test Set Pipeline (`scripts/02_build_test_set.ipynb`)

This notebook mirrors the preprocessing and feature extraction workflow of the training pipeline but is tailored for rigorous, uncompromised model evaluation:

- **Independent GLORI Negative Set:** Incorporates a separate, un-subsampled negative GLORI dataset to preserve natural, negative-heavy class imbalance.
- **Chromosome Isolation:** Restricts the test set **exclusively to `chr1`, `chr17`, and `chr21`** to prevent any spatial genomic or sequence leakage from the training chromosomes.
- **Multi-Model Test Export:** Exports the corresponding test matrices for Models 1, 2, and 3 alongside coordinate tracking files (`test_ids.csv`).


### 3. Sequence Context Extraction (`scripts/extract_sequence_context.R`)

- **Sequence Context Extraction:** An R script utilizing `Biostrings` to extract genomic 5-mer sequence contexts (\(S_{-2}\) to \(S_{+2}\)) from the reference genome (`BSgenome.Hsapiens.UCSC.hg38`), automatically handling reverse-complementation for negative-strand mappings before one-hot encoding.


