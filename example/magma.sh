#!/bin/bash

# Set up directories
export MAGMA_DIR=/mnt/isilon/gandal_lab/mayl/01_GWAS_tools/MAGMA
export DATA=$(pwd)
export OUTPUT=$(pwd)

# Input trait as a parameter
trait=$1  # Specify the trait as the first argument when running the script
sample_size=$2


# Create directories if they do not exist
mkdir -p $DATA
mkdir -p $OUTPUT

# magma annotation
$MAGMA_DIR/magma \
    --snp-loc $DATA/$trait.location \
    --annotate window=5,5 --gene-loc $MAGMA_DIR/NCBI37.3.gene.loc \
    --out $OUTPUT/${trait}_processed.annotation

# gene-based association analysis
$MAGMA_DIR/magma \
    --bfile $MAGMA_DIR/g1000_eur/g1000_eur \
    --pval $DATA/$trait.Pval \
    N=$sample_size \
    --gene-annot $OUTPUT/${trait}_processed.annotation.genes.annot \
    --out $OUTPUT/${trait}_processed_magma_results


