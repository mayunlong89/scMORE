### GWAS Summary Statistics Processing for MAGMA

```bash
# Step 1: Perform quality control on all SNPs based on MAF and sex chromosomes
# In our analysis, SNPs with minor allele frequencies (MAF) < 0.01 or located on sex chromosomes (ChrX and ChrY) were excluded.
#' @param MAF >= 0.01
#' @param Exclude chromosomes ChrX and ChrY

# Step 1: Generate mono_count.location file
# Header: SNP, CHR, POS
cat ieu-b-31.vcf | grep -v "##" | awk '{print $3,$1,$2}' > mono_count.location

# Manually edit the header (e.g., add column names)
vim mono_count.location

# Step 2: Generate mono_count.Pval file

## --- Contains 6 columns: ES:SE:LP:AF:SS:ID
# LP-based P-value (log-transformed P)
#cat ieu-b-31.vcf | grep -v "#" | cut -f 10 | sed 's/:/\t/g' | awk '$4>=0.01' | awk '{print $6,$3}' > mono_count.Pval

# Raw P-value (anti-log-transformation from LP)
cat ieu-b-31.vcf | grep -v "#" | cut -f 10 | sed 's/:/\t/g' | awk '$4>=0.01' | awk '{print $6,10^(-$3)}' > mono_count.Pval

## --- Contains 5 columns: ES:SE:LP:AF:ID
# LP-based P-value (log-transformed P)
#cat ebi-a-GCST004627.vcf | grep -v "#" | cut -f 10 | sed 's/:/\t/g' | awk '$4>=0.01' | awk '{print $5,$3}' > lymp_count.Pval

# Raw P-value (anti-log-transformation from LP)
cat ebi-a-GCST004627.vcf | grep -v "#" | cut -f 10 | sed 's/:/\t/g' | awk '$4>=0.01' | awk '{print $5,10^(-$3)}' > lymp_count.Pval

```

### Example for HDLC_chr3 SNP positions and raw P-values
```bash
zcat ../08_LDSC/BBJ_HDLC.txt.gz | awk 'NR>1 && $2==3 {print $1,$2,$3}' > HDLC_chr3.magma.input.snp.chr.pos.txt
zcat ../08_LDSC/BBJ_HDLC.txt.gz | awk 'NR>1 && $2==3 {print $1,10^(-$11)}' > HDLC_chr3.magma.input.p.txt

# Step 3: Using R to prepare data format for MAGMA
# Convert LP (log-transformed P-value) to raw P-value using anti-log-transformation
```

```R
###----------------------------------------------------------------------------------
library(dplyr)
library(tidyr)

# Read mono_count.Pval file with the header: SNP LP
ids <- read.table("mono_count.Pval", header=F, sep = " ")

# Perform anti-log transformation to convert LP values to raw P-values
ids$Pval <- 10^(-as.numeric(ids$V2))

# Retain only the SNP and P columns
ids_2 <- ids[,c(1,3)]
colnames(ids_2) <- c("SNP","P")

# Write the processed file back to mono_count.Pval
write.table(ids_2, file="mono_count.Pval", row.names=FALSE, col.names=TRUE, quote=FALSE)
###----------------------------------------------------------------------------------
```

```bash
# Filter the processed file to include only SNPs with valid rsIDs
cat mono_count.Pval | grep "^rs" > temp
mv temp mono_count.Pval

# Step 4: Run MAGMA

# Run MAGMA analysis for blood cell traits
srun --pty --mem=15G --time=3-00:00:00 bash magma_bloodcell_count.sh

# Example: Running MAGMA for neutrophil count
srun --pty --mem=15G --time=3-00:00:00 bash ../magma.sh neutr_count 563946

# Note: see 'magma.sh' code in (magma)[https://github.com/mayunlong89/scMORE/blob/main/example/magma.sh]

# Alternatively, run MAGMA directly
bash magma.sh neutr_count 563946
# neutr_count: trait name
# 563946: sample size
```

