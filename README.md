# tcga-sample-classification-and-ttest-analysis
R functions for TCGA RNA-seq sample classification and gene-level statistical testing

This repository contains simple reusable R functions for:

- Assigning tumor vs normal samples from TCGA-style sample IDs
- Performing gene-level t-tests between conditions

## Functions

### 1. get_subtype()
Classifies samples based on TCGA barcode suffix:
- 01, 02, 06 → cancer
- 11 → normal

### 2. run_gene_ttest()
Performs per-gene t-tests between cancer and normal expression datasets.
