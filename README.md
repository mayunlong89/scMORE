
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

#load R packages
library(scMORE)
library(GenomicRanges)
library(IRanges)
library(Seurat)
library(Signac)

#load single-cell data
single_cell <- readRDS("10X_PBMC_downsample_2000cells.rds")

#load GWAS summary data
snp_info <- read.csv("lymphocyte_count_maf0.01.txt",header=T,sep = "\t",stringsAsFactors = FALSE)

#load gene-level association results from MAGMA or FUMA
gene_info <- read.table("lymp_count_processed_magma_results.genes.out",header = TRUE)

#run scMORE 
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



