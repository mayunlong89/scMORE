

##---包含6列： ES:SE:LP:AF:SS:ID
cat ieu-b-31.vcf | grep -v "#" | cut -f 1,2,10 | sed 's/:/\t/g' | awk '$6>=0.01' > monocyte_count_maf0.01.txt
 

#----包含5列：ES:SE:LP:AF:ID
cat ebi-a-GCST004627.vcf | grep -v "#" | cut -f 1,2,10 | sed 's/:/\t/g' | awk '$6>=0.01' > lymphocyte_count_maf0.01.txt

##NOTE: 这里5和6列代码一样,只是提取的Header不一样



##---AD GWAS summary statistics 没有SNP MAF 信息，我利用同样是European Ancestry 样本量更大的PD GWAS summary statistics SNP MAF信息的作为参考信息：去除MAF < 0.01 SNPs.


awk 'NR==FNR {maf[$1]=$2; next} ($1 in maf) && maf[$1]>=0.01' ../02_PD/PD_gwas_maf.txt AD_maf.txt > AD_filtered_maf0.01.txt




