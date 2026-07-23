
<!-- README.md is generated from README.Rmd. Please edit that file -->

<table>
  <tr>
    <td>
      <img src="https://github.com/mayunlong89/scMORE/blob/main/data/scMORE_logo.png" width="250">
    </td>
    <td>
      <h1>scMORE</h1>
      <p>A scalable framework for inferring trait-relevant TF-eRegulons by integrating single-cell multiomic data with GWAS data</p>
    </td>
  </tr>
</table>

<!-- badges: start -->
<!-- badges: end -->

[scMORE](https://github.com/mayunlong89/scMORE) (single cell MultiOmics Regulon Enrichment) is designed to identify disease-relevant eRegulons by leveraging GWAS summary statistics and multimodal single-cell measurements, where gene expression and chromatin accessibility are either measured for each individual cell or integrated into metacell or clusters to capture both modalities.


![Workflow](https://github.com/mayunlong89/scMORE/blob/main/data/Figure%201.png)

## Installation

We recommend installing scMORE via Github using devtools:

``` r
# install packages
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

devtools::install_github("mayunlong89/scMORE")


```
See the DESCRIPTION file for a complete list of R dependencies. If the R dependencies are already installed, installation should finish in a few minutes. You can find the old version from
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
single_cell <- readRDS("10X_PBMC_cells.rds")

# Assign predefined cell types
Idents(single_cell) <- single_cell$cell_type


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
       peak2gene_method = 'Signac',
       infer_method = 'glm',
       method = 'cosine',
       nSeed = 1234)


```



| Functions                           | Description                                                                                                                 |
|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| `scMore()`                         | The main function to fit the model for identifying cell type-specific eRegulons relevant to traits       |
| `createRegulon()`                  | Construct global TF-gene regulotory network using multimodal single-cell measurements                    |
| `getSpecificity()`                 | Calculate the cell type-level specificity score of each gene or TF using the cosine similarity algorithm                    |
| `getGeneScore()`                   | Extract and format gene-based genetic association results from MAGMA or FUMA                                 |
| `regulon2disease()`                | Identify regulons relevant to disease by integrating GWAS-derived gene relevance scores (GRS) and cell type-specificity scores (CTS)     |
| `peak2gene()`                      | Link genomic peaks to genes |
| `snp2peak()`                       | Link SNPs to genomic peaks  |
| `getPeakScore()`                   | Calculate the gene relevance scores (GRS) using the product of peak relevance, peak-gene strength, and gene relevance |
| `getRiskScore()`                   | Extract the highest peak-gene strength scores for each target (node) |
| `getRegulonScore()`                | Calculate trait-associated regulon score (TRS) for each eRegulon |
| `getRandomScore()`                 | Generate random eRegulon scores with matched regulon size for each REAL eRegulon |
| `max_min_scale()`                  | Scale the specificity score between 0 and 1 across all cells                                                                |
| `getEnergyScore()`                 | Calculate the energy score for each dataset analysis (using for method assessment). One can use raw E-score, or use log2(rawEscore+1)+1 to scale the E-score  |




## Assigning cell types to single-cell data

```r

# Assign cell types
Idents(single_cell) <- single_cell$cell_type

```


### Example input format

```r
#1) Single-cell data

The input format of single-cell data: Seurat-generated S4 object based on the GRCh38 genome assembly.

scMORE is fully compatiable with Seurat, a widely-used single-cell analysis tool.

#2) snp_info (GWAS summary statistics)
#CHR: chromosome;
#POS: position;
#ES: effect size;
#SE: standard error;
#LP: -log10(p-value);
#AF: allele frequency;
#SZ: sample size;
#SNP: SNP ID.

CHR	POS	ES	SE	LP	AF	SZ	SNP
1	49298	0.003358	0.011295	0.115634	0.823742	35253	rs10399793
1	54712	0.045321	0.044239	0.514848	0.010379	35253	rs573184866
1	55326	-0.010242	0.033275	0.120193	0.019596	35253	rs3107975


#3) gene_info (MAGMA results)
#GENE: Entrez Gene ID (Gene symbols will be automatically annotated);
#CHR: chromosome;
#START: start position;
#STOP: end position;
#NSNPs: number of SNPs annotated;
#N: sample sizes;
#ZSTAT: Z-scores;
#P: raw P values.

GENE       CHR      START       STOP  NSNPS  NPARAM       N        ZSTAT            P
148398       1     854993     884961     76      20  482730       0.7726      0.21988
26155        1     874583     899679     58      13  482730       0.4058      0.34244
339451       1     890967     906099     34       8  482730      0.70319      0.24097
84069        1     896872     915488     47      16  482730     -0.17594      0.56983

#4) gene_info (FUMA results)
#GENE: ENSEMBL ID;
#CHR: chromosome;
#START: start position;
#STOP: end position;
#NSNPs: number of SNPs annotated;
#N: sample sizes;
#ZSTAT: Z-scores;
#P: raw P values;
#SYMBOL: Gene symbols.

GENE	CHR	START	STOP	NSNPS	NPARAM	N	ZSTAT	P	SYMBOL
ENSG00000237683	1	134901	139379	6	4	171643	0.12078	0.45193	AL627309.1
ENSG00000269831	1	738532	739137	4	2	171643	0.52502	0.29979	AL669831.1
ENSG00000269308	1	818043	819983	6	4	171643	-1.23	0.89065	AL645608.2
ENSG00000187634	1	860260	879955	211	61	171643	-0.61912	0.73208	SAMD11
ENSG00000268179	1	861264	866445	77	32	171643	-0.26316	0.60379	AL645608.1
ENSG00000188976	1	879584	894689	134	37	171643	-0.2079	0.58235	NOC2L
ENSG00000187961	1	895967	901095	66	21	171643	0.24218	0.40432	KLHL17


````
> All GWAS variant coordinates were harmonized to GRCh38 using [UCSC liftover](https://genome.ucsc.edu/cgi-bin/hgLiftOver), and only uniquely mapped variants were retained.
> scATAC-seq peaks, gene annotations, conserved elements and reference sequences were represented in GRCh38 before SNP-to-peak mapping.

### Generate MAGMA results

```shell
# 1) MAGMA codes for generating disease-relevant association scores for TFs and target genes

#DIRECTORY
export MAGMA_DIR=/share/pub/mayl/MAGMA
export DATA=/share/pub/mayl/MAGMA_test
export OUTPUT=/share/pub/mayl/MAGMA_test

#MAGMA annotation:
# By default, a 10 kb window centered on the TSS of a gene is used.

$MAGMA_DIR/magma \
    --snp-loc  $DATA/GWAS_UKBiobank_summary_final.hg19.location  \
    --annotate window=10,10 --gene-loc $MAGMA_DIR/NCBI37.3.gene.loc \
    --out $OUTPUT/GWAS_UKBiobank_summary_final_SNP_Gene_annotation  

# gene-based association analysis:
$MAGMA_DIR/magma \
    --bfile $MAGMA_DIR/1000G_data/g1000_eur \
    --pval $DATA/GWAS_UKBiobank_summary_final.results_Pval \
    N=13239 \
    --gene-annot   $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_annotation.genes.annot  \
    --out $OUTPUT/GWAS_UKBiobank_summary_final_SNP_Gene_Analysis_P


# 2) Processing MAGMA-results: 'magma.genes.out'
# gene_info
magma_results <- read.table("magma.genes.out",header = TRUE)



```
###### For more detailed codes on MAGMA tool, please refer to [GWASTutorial](https://cloufield.github.io/GWASTutorial/09_Gene_based_analysis/), and download the MAGMA tool from [CNCR](https://cncr.nl/research/magma/).
> #### FUMA/MAGMA gene_info: may be based on either hg19 or GRCh38 because it provides gene-level association statistics; coordinate consistency is not essential for this input.


# Install other dependent packages
```r

# If your system did not install some dependent packages: COSG, Pando, ArchR, please install first
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install Seurat from CRAN 
install.packages('Seurat')

# or Install Seurat v5 from Github
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)

# Install ArchR package
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

# Install Pando package
devtools::install_github('quadbio/Pando')

# Install COSG package
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github(repo = 'genecell/COSGR')


```


### Citations
Ma et al., Integrating polygenic signals and single-cell multiomics identifies cell-type-specific regulomes critical for immune- and aging-related diseases, [Nature Aging](https://www.nature.com/articles/s43587-025-01027-5), 2026.

### scMORE reproduce
> For more details on the repeat codes and figure generation scripts, please refer to [scMORE_reproduce](https://github.com/mayunlong89/scMORE_reproduce).
 

## scHOB Database
Human organoids are advanced three-dimensional structures that accurately recapitulate key characteristics of human organ development and functions. Unlike two-dimensional cultures lacking critical cell-cell communications, organoids provides a powerful model for recovering complex cellular dynamics involved in developmental and homeostatic processes. Organoids also allow genetic and pharmacological manipulation in a more physiologically relevant context compared to animal models. Although single-cell sequencing advancements have accelerated their biological and therapeutic use, there has been no systematic platform for unified processing and analysis of organoid-based single-cell multiomics data. 

We thus established scHOB ([single-cell Human Organoid Bank](https://schob.su-lab.org/)), a multi-omic single-cell database, consisting of both scRNA-seq and scATAC-seq data on 10 types of widely-adopted human organoids (i.e., brain, lung, heart, eye, liver & bile duct, pancreas, intestine, kidney, and skin) spanning more than 1.5 million cells with 67 main cell types in 385 samples across 83 distinct protocols. see [Github code](https://github.com/mayunlong89/scHOB/tree/main); see [scHOB Website](https://schob.su-lab.org/).


### Application examples:
1. Ma et al., Integration of human organoids single-cell transcriptomic profiles and human genetics repurposes critical cell type-specific drug targets for severe COVID-19. [Cell Proliferation](https://onlinelibrary.wiley.com/doi/full/10.1111/cpr.13558),2024, and see related [Github codes](https://github.com/mayunlong89/scHuman_organoids_COVID19).
2. Ma et al., Sytematic dissection of pleiotropic loci and critical regulons in excitatory neurons and microglia relevant to neuropsychiatric and ocular diseases, [Translational Psychiatry](https://doi.org/10.1038/s41398-025-03243-4), 2025. preprint version see: [Research Square](https://www.researchsquare.com/article/rs-4514542/v1), 2025.
3. Ma, Y. & Su, J. Linking single-cell multiomics with GWAS to reveal key regulators of disease risk. [Nature Aging](https://doi.org/10.1038/s43587-025-01047-1) 1–2 (2026).



### Other references:
1. [scPagwas](https://www.cell.com/cell-genomics/pdf/S2666-979X(23)00180-5.pdf)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8137370.svg)](https://doi.org/10.5281/zenodo.8137370)

2. [scPagwas_py](https://github.com/mayunlong89/scPagwas_py)
We are preparing to release an upgraded version of scPagwas, named scPagwas2, which introduces enhanced methods for calculating genetically associated genes. By incorporating extensive random calculations, this version offers improved result stability. Additionally, we have addressed issues with synchronizing results across single-cell data and cell-type data. Please note that scPagwas2 requires the use of the scPagwas_main2 function to replace the original scPagwas_main. To accommodate the substantial memory demands of single-cell data calculations in R, we have developed a Python version, scPagwas_py, fully synchronized with scPagwas2.0. We will continue to provide updates to further enhance computational efficiency.

4. Ma et al.,Polygenic regression uncovers trait-relevant cellular contexts through pathway activation transformation of single-cell RNA sequencing data. see [Cell Genomics, 2023](https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00180-5), and see related [Github codes](https://github.com/mayunlong89/scPagwas_main)
5. Li, Ma et al. Integrating microbial GWAS and single-cell transcriptomics reveals associations between host cell populations and the gut microbiome. see [Nature Microbiology, 2025](https://www.nature.com/articles/s41564-025-01978-w), and see related [scBPS](https://zenodo.org/records/15073160) codes: [Zenodo](https://zenodo.org/records/15073160).

