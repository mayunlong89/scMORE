
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scMORE

<!-- badges: start -->
<!-- badges: end -->

scMORE (single cell MultiOmics Regulon Enrichment) is designed to identify disease-relevant regulons by leveraging GWAS summary statistics and multimodal single-cell measurements, where gene expression and chromatin accessibility are either measured for each individual cell or integrated into metacell or clusters to capture both modalities.


![Workflow](https://github.com/mayunlong89/scMORE/blob/main/example/Figure%201.png)

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
| `scMore()`                         | The main function to fit the model for identifying cell type-specific regulons relevant to traits       |
| `createRegulon()`                  | Construct global TF-gene regulotory network using multimodal single-cell measurements                    |
| `getSpecificity()`                 | Calculate the cell type-level specificity score of each gene or TF using the cosine similarity algorithm                    |
| `getGeneScore()`                   | Extract and format gene-based genetic association results from MAGMA or FUMA                                 |
| `regulon2disease()`                | Identify regulons relevant to disease by integrating GWAS-derived gene relevance scores (GRS) and cell type-specificity scores (CTS)     |
| `peak2gene()`                      | Link genomic peaks to genes |
| `snp2peak()`                       | Link SNPs to genomic peaks  |
| `getPeakScore()`                   | Calculate the gene relevance scores (GRS) using the product of peak relevance, peak-gene strength, and gene relevance |
| `getRiskScore()`                   | Extract the highest peak-gene strength scores for each target (node) |
| `getRegulonScore()`                | Calculate trait-associated regulon score (TRS) for each regulon |
| `getRandomScore()`                 | Generate random regulon scores with matched regulon size for each actual regulon |
| `max_min_scale()`                  | Scale the specificity score between 0 and 1 across all cells                                                                |
| `getEnergyScore()`                 | Calculate the energy score for each dataset analysis (using for method assessment) |




## Assigning cell types to single-cell data

```r

# Assign cell types
Idents(single_cell) <- single_cell$cell_type

```


### Example input format

```r
#1) Single-cell data

The input format of single-cell data: Seurat-generated S4 object.

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
1	72526	0.009192	0.023307	0.159069	0.038342	35253	rs547237130
1	79033	0.030693	0.046669	0.291776	0.998702	403994	rs2462495
1	84139	0.026065	0.03179	0.38482	0.020789	35253	rs183605470
1	86028	-0.006518	0.018266	0.141932	0.061284	35253	rs114608975
1	86192	0.00243	0.024492	0.0357588	0.034608	35253	rs548281277
1	91536	0.006541	0.008807	0.339441	0.562561	35253	rs1251109649
1	526736	-0.013042	0.024906	0.221458	0.031698	35253	rs28863004
1	532929	0.029139	0.04165	0.314995	0.012157	35253	rs12240002
1	546697	-0.000347	0.005137	0.0240662	0.913471	403994	rs12025928
1	564862	-0.04748	0.028287	1.03037	0.026161	35253	rs1988726
1	565130	-0.006163	0.042891	0.0527045	0.011434	35253	rs371431021
1	565196	-0.038904	0.041193	0.462261	0.012116	35253	rs139723294
1	565265	-0.008445	0.039132	0.0813753	0.013217	35253	rs868183254
1	565282	0.039097	0.036339	0.549828	0.014143	35253	rs567227003
1	565490	0.005166	0.019022	0.104616	0.059681	35253	rs7349153
1	566024	-0.013942	0.032681	0.174126	0.020699	35253	rs6421779

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
ENSG00000187583	1	901877	911245	103	47	171643	1.4337	0.075832	PLEKHN1
ENSG00000187642	1	910579	917497	48	17	171643	0.87895	0.18971	C1orf170
ENSG00000188290	1	934342	935552	9	5	171643	0.46333	0.32156	HES4
ENSG00000187608	1	948803	949920	10	5	171643	0.086699	0.46546	ISG15
ENSG00000188157	1	955503	991496	337	56	171643	-2.0164	0.97812	AGRN
ENSG00000237330	1	1006346	1009687	28	17	171643	0.10422	0.4585	RNF223
ENSG00000131591	1	1017198	1051741	293	51	171643	-0.65625	0.74417	C1orf159
ENSG00000162571	1	1109264	1133315	286	41	171643	-0.097567	0.53886	TTLL10
ENSG00000186891	1	1138888	1142071	22	9	171643	1.2827	0.099805	TNFRSF18
ENSG00000186827	1	1146706	1149518	19	10	171643	-0.19467	0.57717	TNFRSF4
ENSG00000078808	1	1152288	1167411	182	19	171643	0.4813	0.31515	SDF4
ENSG00000176022	1	1167629	1170421	19	9	171643	1.9922	0.023177	B3GALT6
ENSG00000184163	1	1177826	1182102	28	11	171643	0.19555	0.42248	FAM132A


````



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
    --out $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_annotation  

# gene-based association analysis:
$MAGMA_DIR/magma \
    --bfile $MAGMA_DIR/1000G_data/g1000_eur \
    --pval $DATA/GWAS_UKBiobank_summary_final.results_Pval \
    N=13239 \
    --gene-annot   $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_annotation.genes.annot  \
    --out $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_Analysis_P


# 2) Processing MAGMA-results: 'magma.genes.out'
# gene_info
magma_results <- read.table("magma.genes.out",header = TRUE)



```
#For more detailed codes on MAGMA tool, please refer to [here](https://cloufield.github.io/GWASTutorial/09_Gene_based_analysis/)


# Install other dependent packages
```r

# If your system did not install some dependent packages: COSG, Pando, ArchR, please install first
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install ArchR package
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

# Install Pando package
devtools::install_github('quadbio/Pando')

# Install COSG package
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github(repo = 'genecell/COSGR')


```


### Citations
Ma et al., scMORE uncovers cellular context-specific regulomes linked to complex diseases by integrating polygenic signals and single-cell multiomics, `In submission` 2025



## scHOB Database
Human organoids are advanced three-dimensional structures that accurately recapitulate key characteristics of human organ development and functions. Unlike two-dimensional cultures lacking critical cell-cell communications, organoids provides a powerful model for recovering complex cellular dynamics involved in developmental and homeostatic processes. Organoids also allow genetic and pharmacological manipulation in a more physiologically relevant context compared to animal models. Although single-cell sequencing advancements have accelerated their biological and therapeutic use, there has been no systematic platform for unified processing and analysis of organoid-based single-cell multiomics data. 

We thus established scHOB (single-cell Human Organoid Bank), a multi-omic single-cell database, consisting of both scRNA-seq and scATAC-seq data on 10 types of widely-adopted human organoids (i.e., brain, lung, heart, eye, liver & bile duct, pancreas, intestine, kidney, and skin) spanning more than 1.5 million cells with 67 main cell types in 385 samples across 83 distinct protocols. see [Github code](https://github.com/mayunlong89/scHOB/tree/main); see [scHOB Website](https://schob.su-lab.org/).


### Application examples:
1. Ma et al., Integration of human organoids single-cell transcriptomic profiles and human genetics repurposes critical cell type-specific drug targets for severe COVID-19. [Cell Proliferation](https://onlinelibrary.wiley.com/doi/full/10.1111/cpr.13558),2024, and see related [Github codes](https://github.com/mayunlong89/scHuman_organoids_COVID19).
2. Ma et al., Sytematic dissection of pleiotropic loci and critical regulons in excitatory neurons and microglia relevant to neuropsychiatric and ocular diseases, [Translational Psychiatry](https://doi.org/10.1038/s41398-025-03243-4), 2025. preprint version see: [Research Square](https://www.researchsquare.com/article/rs-4514542/v1), 2024.


### Other references:
1. [scPagwas](https://www.cell.com/cell-genomics/pdf/S2666-979X(23)00180-5.pdf)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8137370.svg)](https://doi.org/10.5281/zenodo.8137370)

2. [scPagwas_py](https://github.com/mayunlong89/scPagwas_py)
We are preparing to release an upgraded version of scPagwas, named scPagwas2, which introduces enhanced methods for calculating genetically associated genes. By incorporating extensive random calculations, this version offers improved result stability. Additionally, we have addressed issues with synchronizing results across single-cell data and cell-type data. Please note that scPagwas2 requires the use of the scPagwas_main2 function to replace the original scPagwas_main. To accommodate the substantial memory demands of single-cell data calculations in R, we have developed a Python version, scPagwas_py, fully synchronized with scPagwas2.0. We will continue to provide updates to further enhance computational efficiency.

4. Ma et al.,Polygenic regression uncovers trait-relevant cellular contexts through pathway activation transformation of single-cell RNA sequencing data. see [Cell Genomics, 2023](https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00180-5), and see related [Github codes](https://github.com/mayunlong89/scPagwas_main)

