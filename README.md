
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scMORE

<!-- badges: start -->
<!-- badges: end -->

scMORE is designed to identify disease-relevant regulons by leveraging GWAS summary statistics and multimodal single-cell measurements, where gene expression and chromatin accessibility are either measured for each individual cell or integrated into metacell or clusters to capture both modalities.

## Installation

We recommend installing scMORE via Github using devtools:

``` r
# install.packages("devtools")
devtools::install_github("mayunlong89/scMORE")
```
See the DESCRIPTION file for a complete list of R dependencies. If the R dependencies are already installed, installation should finish in a few minutes. You can install the old version of scMORE from
[GitHub](https://github.com/mayunlong89/ctDRTF).


## How to run scMORE
```r

library(scMORE)

#' @title Single-Cell MultiOmics REgulon discovery (scMORE)
#'
#' @description
#' Main function for identifying disease-relevant regulons using single-cell data.
#'
#' @param single_cell The input single-cell data, typically a Seurat or SingleCellExperiment object.
#' @param snp_info GWAS summary statistics. Must include the following columns:
#'                 `chr` (chromosome), `pos` (position), `rsID` (SNP ID), `beta` (effect size),
#'                 and `p-value` (significance level).
#' @param n_targets Minimum number of target genes required in a given regulon. Default: 5.
#' @param gene_info Gene-based genetic association results from MAGMA or FUMA.
#'                  The function supports results in standard MAGMA/FUMA formats.
#' @param perm_n Number of Monte Carlo permutations to perform for significance testing. Default: 1000.
#' @param theta Weighting factor to adjust the specificity score of target genes within each regulon.
#'              Range: 0.1 to 1. Default: 0.5.
#' @param alpha Flexibility parameter for penalization in the inference model. Default: 1.
#' @param buffer Numeric value specifying the flanking region size for genomic peaks.
#'               Extends the peak range upstream and downstream. For example, `buffer = 500`
#'               extends 500 bp on both sides of a peak. Default: 500 bp.
#' @param p1 Threshold for statistical significance of the cell type-specificity score
#'           for each regulon. Default: 0.05.
#' @param p2 Threshold for statistical significance of the genetic risk score
#'           for each regulon. Default: 0.05.
#' @param p3 Threshold for statistical significance of the trait-associated regulon score (TARS).
#'           Default: 0.05.
#' @param peak2gene_method Method for mapping peaks to genes. Options: `'Signac'` or `'GREAT'`.
#' @param infer_method Method for inference modeling. Options:
#'                     - Generalized Linear Models (`'glm'`)
#'                     - Regularized GLMs (`'glmnet'`, `'cv.glmnet'`)
#'                     - Bayesian regression models (`'brms'`).
#' @param top_n Number of top targets for each TF used to calculate its importance. Default: 5.
#' @param nSeed Random seed for reproducibility (i.e., `set.seed()`). Default: 1234.
#'
#' @return A list containing:
#'         - Disease-relevant regulons
#'         - Their genetic and cell-type associations
#'
#' @export
#'

scMore(single_cell,
       snp_info,
       gene_info,
       n_targets = 5,
       perm_n = 1000,
       theta = 0.5,
       alpha = 1,
       buffer = 500,
       top_n = 5,
       p1 = 0.05,
       p2 = 0.05,
       p3 = 0.05,
       nSeed = 1234)



```


| Function                           | Description                                                                                                                 |
|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| `scMore()`                         | The main function to fit the model for identifying cell type-specific regulons relevant to traits       |
| `createRegulon()`                   Construct global TF-gene regulotory network using multimodal single-cell measurements                    |
| `COSR_pre_func()`                  | Calculate the cell type-level specificity score of each gene or TF using the cosine similarity algorithm                    |
| `COSR_func_weight()`               | Identify cell type-specific regulons relevant to disease using polygenic enrichment method                                  |
| `MC_JSI_score_func_weight()`       | Calculate the empirical P value for each regulon-disease association using Monote Carlo permutation algorithm               |
| `max_min_scale()`                  | Scale the specificity score between 0 and 1 across all cells                                                                |



