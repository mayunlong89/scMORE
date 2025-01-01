
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scMORE

<!-- badges: start -->
<!-- badges: end -->

scMORE is designed to identify disease-relevant regulons by leveraging GWAS summary statistics and multimodal single-cell measurements, where gene expression and chromatin accessibility are either measured for each individual cell or integrated into metacell or clusters to capture both modalities.


![Workflow](https://github.com/mayunlong89/scMORE/blob/main/example/Figure%201.png)

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


| Functions                           | Description                                                                                                                 |
|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| `scMore()`                         | The main function to fit the model for identifying cell type-specific regulons relevant to traits       |
| `createRegulon()`                  | Construct global TF-gene regulotory network using multimodal single-cell measurements                    |
| `getSpecificity()`                 | Calculate the cell type-level specificity score of each gene or TF using the cosine similarity algorithm                    |
| `getGeneScore()`                   | Extract and format gene-based genetic association results from MAGMA or FUMA                                 |
| `regulon2disease()`                | Identify regulons relevant to disease by integrating GWAS-derived gene relevance scores (GRS) and cell type-specificity scores (CTS)     |
| `peak2gene()`                      | Link genes to genomic peaks |
| `getPeakScore()`                   | Calculate the gene relevance scores (GRS) using the product of peak relevance, peak-gene strength, and gene relevance |
| `getRiskScore()`                   | Extract the highest peak-gene strength scores for each target (node) |
| `getRegulonScore()`                | Calculate trait-associated regulon score (TRS) for each regulon |
| `getRandomScore()`                 | Generate random regulon scores with matched regulon size for each actual regulon |
| `max_min_scale()`                  | Scale the specificity score between 0 and 1 across all cells                                                                |



## Assigning cell types to single-cell data

```r

Idents(single_cell) <- single_cell$cell_type

```

### Example format

```
1) Single-cell data

The input format of single-cell data: Seurat-generated S4 object.

scMORE is fully compatiable with Seurat, a widely-used single-cell analysis tool.


2) gene_info (MAGMA results)
GENE       CHR      START       STOP  NSNPS  NPARAM       N        ZSTAT            P
148398       1     854993     884961     76      20  482730       0.7726      0.21988
26155        1     874583     899679     58      13  482730       0.4058      0.34244
339451       1     890967     906099     34       8  482730      0.70319      0.24097
84069        1     896872     915488     47      16  482730     -0.17594      0.56983
84808        1     905579     922473     56      12  482730    -0.077128      0.53074
57801        1     929342     941608     46       9  482730      0.40403      0.34309
9636         1     943847     954920     28       6  482730       1.2899     0.098549
375790       1     950503     996499    112      17  482730       0.1929      0.42352
401934       1    1002126    1014687     19       4  482730     -0.76977      0.77928
54991        1    1012198    1056736    119      14  482730       1.0179      0.15436
254173       1    1104286    1138315    154      22  482730     -0.84448       0.8008
8784         1    1133888    1147163     38       7  482730     0.098027      0.46096
7293         1    1141706    1154703     45       8  482730      0.41393      0.33946
51150        1    1147288    1172447    115      10  482730      0.37934      0.35222
126792       1    1162629    1175421     50       8  482730      0.43793      0.33072
388581       1    1172826    1187102     44       8  482730      0.38135      0.35147
118424       1    1184292    1214234     82      14  482730      0.30781      0.37911
6339         1    1210816    1232409     69      19  482730       -1.267      0.89742

3) gene_info (FUMA results)


````


### Generate MAGMA results

```shell
1) MAGMA codes for generating disease-relevant genes

#DIRECTORY
export MAGMA_DIR=/share/pub/mayl/MAGMA
export DATA=/share/pub/mayl/MAGMA_test
export OUTPUT=/share/pub/mayl/MAGMA_test

#MAGMA annotation:

$MAGMA_DIR/magma \
    --snp-loc  $DATA/GWAS_UKBiobank_summary_final.hg19.location  \
    --annotate window=20,20 --gene-loc $MAGMA_DIR/NCBI37.3.gene.loc \
    --out $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_annotation  

#gene-based association analysi:
$MAGMA_DIR/magma \
    --bfile $MAGMA_DIR/1000G_data/g1000_eur \
    --pval $DATA/GWAS_UKBiobank_summary_final.results_Pval \
    N=13239 \
    --gene-annot   $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_annotation.genes.annot  \
    --out $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_Analysis_P


2) Processing MAGMA-results: 'magma.genes.out'
#MAGMA_GWAS_data: all MAGMA-based associations results ranked by -log10(P)
#header of MAGMA_GWAS_data: SYMBOL, logP, ZSTAT

magma_results <- read.table("magma.genes.out",header = TRUE)
magma_results <- magma_results %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
MAGMA_GWAS_data <- magma_results[,c(10,11,8)]


```
#For more detailed codes on MAGMA tool, please refer to [here](https://cloufield.github.io/GWASTutorial/09_Gene_based_analysis/)


# Citations
1. Ma et al., Sytematic dissection of pleiotropic loci and critical regulons in exhibitory neurons and microglia relevant to neuropsychiatric and ocular diseases, Translational Psychiatry, 2025. preprint version see: [Research Square](https://www.researchsquare.com/article/rs-4514542/v1), 2024.
2. Ma et al., Polygenic network enrichment identifies cellular context-specific regulons relevant to diseases by integration of single-cell multiomic data, `Genome Biology`,(under review), 2024


# scHOB Database
Human organoids are advanced three-dimensional structures that accurately recapitulate key characteristics of human organ development and functions. Unlike two-dimensional cultures lacking critical cell-cell communications, organoids provides a powerful model for recovering complex cellular dynamics involved in developmental and homeostatic processes. Organoids also allow genetic and pharmacological manipulation in a more physiologically relevant context compared to animal models. Although single-cell sequencing advancements have accelerated their biological and therapeutic use, there has been no systematic platform for unified processing and analysis of organoid-based single-cell multiomics data. 

We thus established scHOB (single-cell Human Organoid Bank), a multi-omic single-cell database, consisting of both scRNA-seq and scATAC-seq data on 10 types of widely-adopted human organoids (i.e., brain, lung, heart, eye, liver & bile duct, pancreas, intestine, kidney, and skin) spanning more than 1.5 million cells with 67 main cell types in 385 samples across 83 distinct protocols. see [Github code](https://github.com/mayunlong89/scHOB/tree/main); see [scHOB Website](https://schob.su-lab.org/).
The single-cell multiome data in scHOB have been used by ctDRTF, see Ma et al. Translational Psychiatry, 2024.

# Application example of scHOB database:
Ma et al., Integration of human organoids single-cell transcriptomic profiles and human genetics repurposes critical cell type-specific drug targets for severe COVID-19. [Cell Proliferation](https://onlinelibrary.wiley.com/doi/full/10.1111/cpr.13558),2024, and see related [Github codes](https://github.com/mayunlong89/scHuman_organoids_COVID19).


# Other references:
 
1. [scPagwas](https://www.cell.com/cell-genomics/pdf/S2666-979X(23)00180-5.pdf)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8137370.svg)](https://doi.org/10.5281/zenodo.8137370)

2. Development of novel polygenic regression method scPagwas for integrating scRNA-seq data with GWAS on complex diseases. see [Ma et al. Cell Genomics, 2023](https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00180-5), and see related [Github codes](https://github.com/mayunlong89/scPagwas_main)

