
## This is the code for scMORE to link SNPs to peaks using fine-mapped SNPs with high posterior probabilities.

#### 2025-08-06


### Interaction between GWAS Signals and Chromatin Accessibility Peaks
```R
setwd('/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/sigpeaks')
library(data.table)


gwas <- fread('/share2/pub/zhouyj/zhouyj/ctDRTF_test/Method_envalue/GWAS/Lym_cp_zscore/lymphocyte_count_maf0.01.txt', header = TRUE, sep = '\t')
gwas$CHROM <- as.character(gwas$CHROM)

output_dir <- '/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/sigpeaks/peak_snp'
#dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
peak_file <- '/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/sigpeaks/SigGRN.500bp.bed'
peaks <- read.table(peak_file, header = FALSE, sep = '\t')
colnames(peaks) <- c("chr", "start", "end")
peaks$chr <- as.character(peaks$chr)
for (i in 1:nrow(peaks)) {
  peak_chr <- peaks$chr[i]
  peak_start <- peaks$start[i]
  peak_end <- peaks$end[i]


  peak_start_file <- peak_start + 500
  peak_end_file <- peak_end - 500


  snps_in_region <- subset(gwas, CHROM == peak_chr & POS >= peak_start & POS <= peak_end)

  if (nrow(snps_in_region) > 0) {
    file_name <- paste0(output_dir, "/peaks_", peak_chr, "_", peak_start_file, "_", peak_end_file, ".txt")
    write.table(snps_in_region, file = file_name, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}



```
