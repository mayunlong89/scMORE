
#### 1) MAGMA codes for generating disease-relevant association scores for TFs and target genes

```R
#DIRECTORY
export MAGMA_DIR=/share/pub/mayl/MAGMA
export DATA=/share/pub/mayl/MAGMA_test
export OUTPUT=/share/pub/mayl/MAGMA_test

#MAGMA annotation:
#By default, a 10 kb window centered on the TSS of a gene is used.

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

gene_info <- magma_results

```

#### Input format for MAGMA
```R
# Using PD GWAS as an example.
#1) PD.location
SNP CHR POS
rs561234294 1 54591
rs2462492 1 54676
rs534350410 1 79188

#2) PD.Pval (no header: SNP, P-value)
rs561234294 0.012
rs2462492 0.0003
rs534350410 0.096269

```

#### Output format of MAGMA

```R
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
```





