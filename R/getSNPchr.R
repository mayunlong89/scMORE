#'
#' @title Format SNP chromosome names
#'
#' @param snp_bed Data frame of genetic association results in .bed format
#' @param getSNPchr Logical, whether to add "chr" prefix to the chromosome column (default: TRUE)
#' @return Updated data frame with "chr" prefix added to the chromosome column if missing
#' @export
#'
getSNPchr <- function(snp_bed) {
  # Check if the first column contains "chr"
  if (!all(grepl("^chr", snp_bed[[1]]))) {
    # If "chr" prefix is missing, add it
    snp_bed[[1]] <- paste0("chr", snp_bed[[1]])
  }
  return(snp_bed)
}

