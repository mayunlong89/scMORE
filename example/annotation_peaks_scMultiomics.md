
# Annotate Peaks with Gene Names for Single-cell Multi-omics Data

> This script demonstrates how to annotate peak information with gene symbols using GENCODE GTF and biomaRt, and how to integrate it into a Seurat object for downstream regulon analysis.
> 

## Step 1: Load Required Libraries
```r
library(GenomicFeatures)
library(biomaRt)
library(GenomicRanges)
```

## Step 2: Load GENCODE GTF and Create TxDb
```r
gtf_file <- "/Users/mayunlong/Desktop/UPENN/genome_references/gencode.v44.annotation.gtf"

# Create TxDb object from GTF
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# Extract gene annotations
genes <- genes(txdb)
genome(genes) <- "hg38"  # Set genome build
```
> Note: users could download the gencode.v44.annotation.gtf as follow:
```bash
# Download for Unix system

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip gencode.v44.annotation.gtf.gz

```

## Step 3: Annotate Gene IDs with HGNC Symbols
```r
# Remove version numbers from gene IDs (e.g., ENSG000001234.5 -> ENSG000001234)
gene_ids <- gsub("\\..*", "", genes$gene_id)

# Connect to Ensembl biomart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get HGNC gene symbols
annot <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
               filters = "ensembl_gene_id",
               values = gene_ids,
               mart = mart)

# Match annotation to gene order
idx <- match(gene_ids, annot$ensembl_gene_id)
gene_names <- annot$hgnc_symbol[idx]

# Add gene_name to genes metadata
genes$gene_name <- gene_names
```

## Step 4: Assign Annotation to Seurat Object
```r
# Annotate peak ranges with gene info
Annotation(bmmc_subset[["peaks"]]) <- genes
```
>💡 This step attaches gene information (including HGNC symbols) to your peaks, which will be used for downstream analyses.


## Step 5: Run Regulon Inference (e.g., using Pando)
```r
# Set RNA as active assay
DefaultAssay(bmmc_subset) <- "RNA"

# Compute regulons using peak-gene relationships
grn_outputs_bmmc <- createRegulon(bmmc_subset, n_targets = 5)
```
>To run createRegulon() on new single-cell multi-omics data (e.g., RNA + ATAC), make sure that peak-to-gene annotations are properly added. The above pipeline shows how to use GENCODE and Ensembl biomaRt to map peaks to gene symbols and assign them to your Seurat object using Annotation(). This is a required step before downstream regulon or GRN analysis.





