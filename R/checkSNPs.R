#'
#' @title Check and format SNP chromosome names
#'
#' @description
#' This function checks if the required columns are present in a given SNP data frame,
#' and optionally adds a "chr" prefix to the chromosome column if specified.
#'
#' @param snp_info Data frame of genetic association results.
#'        -required_columns A character vector of required column names. Default: c("CHR", "POS", "LP", "SNP").
#' @param add_chr_prefix Logical. Whether to add "chr" prefix to the "CHR" column if not already present. Default: FALSE.
#' @return The modified data frame if all required columns are present. Otherwise, an error is thrown.
#' @export
checkSNPs <- function(snp_info, add_chr_prefix = TRUE) {
  # Ensure snp_info is a data frame
  if (!is.data.frame(snp_info)) {
    stop("Error: The input snp_info must be a data frame.")
  }
  required_columns <- c("CHR", "POS", "LP", "SNP")
  # Check for missing columns
  missing_columns <- setdiff(required_columns, colnames(snp_info))
  if (length(missing_columns) > 0) {
    stop("Error: Missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  # Optionally add "chr" prefix to the CHR column
  if (add_chr_prefix) {
    snp_info$CHR <- ifelse(grepl("^chr", snp_info$CHR, ignore.case = TRUE),
                           as.character(snp_info$CHR),
                           paste0("chr", snp_info$CHR))
  }

  return(snp_info)
}
