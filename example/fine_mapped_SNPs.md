
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

### Constructing 3-Mb Sentinel Regions for Peaks
```R
input_dir <- "/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/sigpeaks"
output_dir <- file.path(input_dir, "sentinel_regions")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(paste0(input_dir, "/peak_snp"), full.names = TRUE)
for (file in files) {
  snps <- read.table(file, header = TRUE, sep = '\t')

  snps <- snps[order(-snps$Zscore), ]

  regions <- data.frame()
  used <- rep(FALSE, nrow(snps))

  for (i in 1:nrow(snps)) {
    if (!used[i]) {
      center <- snps$POS[i]
      chr_label <- snps$CHROM[i]
      start <- center - 1500000
      end <- center + 1500000
      regions <- rbind(regions, data.frame(
        CHROM = chr_label,
        START = start,
        END = end,
        SNP = paste0(chr_label, "_", center)
      ))
      used <- used | (snps$CHROM == chr_label & snps$POS >= start & snps$POS <= end)
    }
  }

  filename <- basename(file)
  output_file <- sub("^peaks_", "sentinel_regions_", filename)
  output_path <- file.path(output_dir, output_file)

  write.table(regions, output_path, row.names = FALSE, quote = FALSE, sep = "\t")
}
```

### Generating LDstore Z-Files from GWAS Summary Statistics
```R
region_dir_base <- output_dir
eur_snp_dir <- "/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/LD_file"  # 每条染色体 SNP list 放在这里

out_dir1 <- "/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/sigpeaks/LD_file/LD_Z_inputs"
out_dir2 <- "/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/sigpeaks/finemap/inputs"

gwas <- gwas[nchar(gwas$A1) == 1 & nchar(gwas$A2) == 1, ]

region_files <- list.files(region_dir_base,  full.names = TRUE)

for (region_file in region_files) {
  regions <- fread(region_file, header = TRUE, sep = "\t")
  regions$CHROM <- as.character(regions$CHROM)

  eur_snp_file <- file.path(eur_snp_dir, paste0("chr", regions$CHROM, ".snp.list"))
  eur_snps <- fread(eur_snp_file, header = TRUE, sep = "\t")
  eur_snps <- eur_snps[nchar(eur_snps$first_allele) == 1 & nchar(eur_snps$alternative_alleles) == 1, ]
  rsid_list <- eur_snps$rsid

  base_name <- tools::file_path_sans_ext(basename(region_file))

  for (i in 1:nrow(regions)) {
    chr <- regions$CHROM[i]
    start <- regions$START[i]
    end <- regions$END[i]

subRegion <- gwas[gwas$CHROM == chr & gwas$POS >= start & gwas$POS <= end, ]

subRegion1 <- data.frame(
  rsid = paste(subRegion$CHROM, subRegion$POS, sep = "_"),
  chromosome = subRegion$CHROM,
  position = subRegion$POS,
  allele1 = subRegion$A1,
  allele2 = subRegion$A2
)
subRegion1 <- subRegion1[subRegion1$rsid %in% rsid_list, ]

subRegion2 <- data.frame(
  rsid = paste(subRegion$CHROM, subRegion$POS, sep = "_"),
  chromosome = subRegion$CHROM,
  position = subRegion$POS,
  allele1 = subRegion$A1,
  allele2 = subRegion$A2,
  maf = ifelse(subRegion$AF <= 0.5, subRegion$AF, 1 - subRegion$AF),
  beta = subRegion$ES,
  se = subRegion$SE
)
subRegion2 <- subRegion2[subRegion2$rsid %in% rsid_list, ]

# === 写出文件 ===
out_file1 <- file.path(out_dir1, paste0(base_name, "_region", i, ".z"))
out_file2 <- file.path(out_dir2, paste0(base_name, "_region", i, ".z"))

write.table(subRegion1, file = out_file1, row.names = FALSE, quote = FALSE, sep = " ")
write.table(subRegion2, file = out_file2, row.names = FALSE, quote = FALSE, sep = " ")

  }
}
```
### ldstore z file 
```bash
region_dir <- "/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/sigpeaks/LD_file/LD_Z_inputs"
output_master <- "/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/sigpeaks/LD_file/LD_masters"
output_ld <- "/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/sigpeaks/LD_file/LD_scores"


#region_files <- list.files(region_dir, pattern = "^sentinel_regions_*\\.z$", full.names = TRUE)
region_files <- list.files(region_dir, full.names = T)


master_lines <- c("z;bgen;bgi;sample;bdose;bcor;ld;n_samples")


for (f in region_files) {

  # 43135550_43137417

  file_base <- sub(".*_(\\d+_\\d+_\\d+)_.*", "\\1", basename(f))
  #"10_60446963_60446999"
  chr <- strsplit(file_base, "_")[[1]][1]
  z <- f
  bgen <- paste0("/share/pub1/guiyy/guiyy/scMORE/vcf_eur/bgen/chr", chr, ".pos.bgen")
  bgi <- paste0(bgen, ".bgi")
  sample <- "/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/LD_file/EUR_vcfs/EUR.data.sample"
  bdose <- paste0(output_ld, "/peaks_", file_base, ".bdose")
  bcor <- paste0(output_ld, "/peaks_", file_base, ".bcor")
  ld <- paste0(output_ld, "/peaks_", file_base, ".ld")
  n_samples <- "503"

  line <- paste(z, bgen, bgi, sample, bdose, bcor, ld, n_samples, sep = ";")
  master_lines <- c(master_lines, line)
  per_master_file <- file.path(output_master, paste0("peaks_", file_base, ".master"))
  writeLines(c("z;bgen;bgi;sample;bdose;bcor;ld;n_samples", line), con = per_master_file)
}

# computing LD matrix

/share2/pub/yaoyh/yaoyh/software/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 --in-files "$master_file" --write-text --read-only-bgen

# Preparing FINEMAP inputs

z_dir <- "/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/sigpeaks/finemap/inputs"
ld_dir <- "/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/sigpeaks/LD_file/LD_scores"
output_master <- "/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/sigpeaks/finemap/master"
output_finemap <- "/share2/pub/yaoyh/yaoyh/1Project/scMORE_revision/FINEMAP/sigpeaks/finemap/outputs"
z_files <- list.files(z_dir, full.names = T)

master_lines <- c("z;ld;snp;config;cred;log;k;n_samples")

for (f in z_files) {
  file_base <- sub(".*_(\\d+_\\d+_\\d+)_.*", "\\1", basename(f))
  #19946918_19948789_region1
  z <- f
  ld <- paste0(ld_dir, "/peaks_", file_base, ".ld")
  snp <- paste0(output_finemap, "/peaks_", file_base, ".snp")
  config <- paste0(output_finemap, "/peaks_", file_base, ".config")
  cred <- paste0(output_finemap, "/peaks_", file_base, ".cred")
  log <- paste0(output_finemap, "/peaks_", file_base, ".log")
  k <- paste0(output_finemap, "/peaks_", file_base, ".k")
  n_samples <- "171643"

  line <- paste(z, ld, snp, config, cred, log, k, n_samples, sep = ";")
  master_lines <- c(master_lines, line)
  per_master_file <- file.path(output_master, paste0("peaks_", file_base, ".finemap.master"))
  writeLines(c("z;ld;snp;config;cred;log;k;n_samples", line), con = per_master_file)
}

# performing fine-mapping

"$FINEMAP_EXEC" --sss --in-files "$master_file" --dataset 1

```








