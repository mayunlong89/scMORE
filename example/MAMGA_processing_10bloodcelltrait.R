
##GWAS Summary processing for MAGMA

1) 首先，所有的SNP MAF and sex chromosome 质控：
#In our applications, SNPs with minor alelle frequencies smaller than 0.01 or on the sex chromosomes (ChrX-Y) were removed.
#' @param MAF >= 0.01
#' @param excluded chrX-Y

#Step 1 for mono_count.location
#Header: SNP, CHR, POS

cat ieu-b-31.vcf | grep -v "##" | awk '{print $3,$1,$2}' > mono_count.location

#Manually modify the header

vim mono_count.location

#Step 2 for mono_count.Pval

##---包含6列： ES:SE:LP:AF:SS:ID
#LP
#cat ieu-b-31.vcf | grep -v "#" | cut -f 10 | sed 's/:/\t/g' | awk '$4>=0.01' | awk '{print $6,$3}' > mono_count.Pval

#Raw P
cat ieu-b-31.vcf | grep -v "#" | cut -f 10 | sed 's/:/\t/g' | awk '$4>=0.01' | awk '{print $6,10^(-$3)}' > mono_count.Pval

#----包含5列：ES:SE:LP:AF:ID
#LP
#cat ebi-a-GCST004627.vcf | grep -v "#" | cut -f 10 | sed 's/:/\t/g' | awk '$4>=0.01' | awk '{print $5,$3}' > lymp_count.Pval

#Raw p
cat ebi-a-GCST004627.vcf | grep -v "#" | cut -f 10 | sed 's/:/\t/g' | awk '$4>=0.01' | awk '{print $5,10^(-$3)}' > lymp_count.Pval


zcat ../08_LDSC/BBJ_HDLC.txt.gz | awk 'NR>1 && $2==3 {print $1,$2,$3}' > HDLC_chr3.magma.input.snp.chr.pos.txt
zcat ../08_LDSC/BBJ_HDLC.txt.gz | awk 'NR>1 && $2==3 {print $1,10^(-$11)}' >  HDLC_chr3.magma.input.p.txt


#Using R to pre-arrange data format for MAGMA 
# anti-log-transformation: --> LP convert to Raw P
>R
###----------------------------------------------------------------------------------
library(dplyr)
library(tidyr)

#mono_count.Pval file header: SNP LP
ids <- read.table("mono_count.Pval", header=F, sep = " ")

#anti-log transformation LPval to raw Pval.
ids$Pval <- 10^(-as.numeric(ids$V2))
                
ids_2 <- ids[,c(1,3)]
colnames(ids_2) <- c("SNP","P")

write.table(ids_2, file="mono_count.Pval", row.names=FALSE, col.names=TRUE, quote=FALSE)

###----------------------------------------------------------------------------------


cat mono_count.Pval | grep "^rs" > temp
mv temp mono_count.Pval
             
                
#run MAGMA

srun --pty --mem=15G --time=3-00:00:00 bash magma_bloodcell_count.sh

srun --pty --mem=15G --time=3-00:00:00 bash  ../magma.sh neutr_count 563946

bash magma.sh neutr_count 563946
# neutr_count: trait name
# 563946: sample size 

