
#### Processing SNPs from GWAS summary statistics to form the 'snp_info' file

```bash
##--- Contains 6 columns: ES:SE:LP:AF:SS:ID
cat ieu-b-31.vcf | grep -v "#" | cut -f 1,2,10 | sed 's/:/\t/g' | awk '$6>=0.01' > monocyte_count_maf0.01.txt

#---- Contains 5 columns: ES:SE:LP:AF:ID
cat ebi-a-GCST004627.vcf | grep -v "#" | cut -f 1,2,10 | sed 's/:/\t/g' | awk '$6>=0.01' > lymphocyte_count_maf0.01.txt

## NOTE: The code for extracting columns 5 and 6 is identical, but the difference lies in the specific column headers being extracted.

##--- AD GWAS summary statistics lack SNP MAF information. To filter SNPs with MAF < 0.01,
##--- We used SNP MAF data from PD GWAS summary statistics (from European ancestry samples with a larger sample size) as a reference.

awk 'NR==FNR {maf[$1]=$2; next} ($1 in maf) && maf[$1]>=0.01' ../02_PD/PD_gwas_maf.txt AD_maf.txt > AD_filtered_maf0.01.txt

## for frailty index GWAS------
#version 1
awk '$6 >= 0.01 {print $2, $3, $7, $8, -log($9)/log(10), $6, $1}' GCST90020053_GRCh37.tsv > frailty_index_maf0.01.txt

#version 2
cat GCST90020053_GRCh37.tsv | awk '$6 >= 0.01' | awk '{print $2, $3, $7, $8, -log($9)/log(10), $6, $1}'  > frailty_index_maf0.01.txt

#change the header into: CHR, POS, ES, SE, LP, AF, SNP


## for lifespan GWAS------
#version 1
awk '$8 >= 0.01 {print $3, $4, $9, $10, -log($11)/log(10), $8, $1}' lifegen_phase2_bothpl_alldr_2017_09_18.tsv > lifespan_maf0.01.txt

#version 2
cat lifegen_phase2_bothpl_alldr_2017_09_18.tsv | awk '$8 >= 0.01' | awk '{print $3, $4, $9, $10, -log($11)/log(10), $8, $1}'  > lifespan_maf0.01.txt

#change the header into: CHR, POS, ES, SE, LP, AF, SNP


## for EAA GWAS------
#version 1
awk '$6 >= 0.01 {print $2, $3, $7, $8, -log($9)/log(10), $6, $1}' GCST90014290_buildGRCh37.tsv > EAA_maf0.01.txt

#version 2
cat GCST90014290_buildGRCh37.tsv | awk '$6 >= 0.01' | awk '{print $2, $3, $7, $8, -log($9)/log(10), $6, $1}'  > EAA_maf0.01.txt

#change the header into: CHR, POS, ES, SE, LP, AF, SNP

## for healthspan GWAS------
#version 1

cat healthspan_summary.csv | sed 's/,/\t/g' | awk '$6 >= 0.01' | awk '{print $2, $3, $7, $8, -log($10)/log(10), $6, $1}'  > healthspan_maf0.01.txt

#change the header into: CHR, POS, ES, SE, LP, AF, SNP

## for mvAge GWAS------
#version 1
awk '$4 >= 0.01 {print $2, $3, $7, $8, -log($9)/log(10), $4, $1}' mvAge.summary.EUR.txt > mvAge_maf0.01.txt

#version 2
cat mvAge.summary.EUR.txt | awk '$4 >= 0.01' | awk '{print $2, $3, $7, $8, -log($9)/log(10), $4, $1}'  > mvAge_maf0.01.txt

#change the header into: CHR, POS, ES, SE, LP, AF, SNP
```
