##--- Contains 6 columns: ES:SE:LP:AF:SS:ID
cat ieu-b-31.vcf | grep -v "#" | cut -f 1,2,10 | sed 's/:/\t/g' | awk '$6>=0.01' > monocyte_count_maf0.01.txt

#---- Contains 5 columns: ES:SE:LP:AF:ID
cat ebi-a-GCST004627.vcf | grep -v "#" | cut -f 1,2,10 | sed 's/:/\t/g' | awk '$6>=0.01' > lymphocyte_count_maf0.01.txt

## NOTE: The code for extracting columns 5 and 6 is identical, but the difference lies in the specific column headers being extracted.

##--- AD GWAS summary statistics lack SNP MAF information. To filter SNPs with MAF < 0.01,
##--- We used SNP MAF data from PD GWAS summary statistics (from European ancestry samples with a larger sample size) as a reference.

awk 'NR==FNR {maf[$1]=$2; next} ($1 in maf) && maf[$1]>=0.01' ../02_PD/PD_gwas_maf.txt AD_maf.txt > AD_filtered_maf0.01.txt
