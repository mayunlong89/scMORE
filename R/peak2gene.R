#'
#' @title Link genes to genomic peaks
#'
#' @param grn_outputs A list containing GRN output data, with a 'Regions' column specifying genomic peak ranges.
#' @return A data frame mapping genes to genomic peaks with associated log-transformed p-values (Strength).
#' @import dplyr
#' @import tidyr
#' @export
#'
peak2gene <- function(grn_outputs) {

  # Extract GRN data
  grn_data <- grn_outputs$grn

  # Ensure the 'Regions' column is of character type
  if (!"Regions" %in% colnames(grn_data)) {
    stop("'Regions' column is missing in the GRN data")
  }
  grn_data$Regions <- as.character(grn_data$Regions)

  # Ensure the 'Pval' column exists for computing Strength
  if (!"Pval" %in% colnames(grn_data)) {
    stop("'Pval' column is missing in the GRN data")
  }

  # Split the 'Regions' column by ";" and expand into multiple rows
  expanded_grn_data <- grn_data %>%
    tidyr::separate_rows(Regions, sep = ";") %>%   # Split rows by ';'
    dplyr::mutate(Strength = -log10(Pval)) %>%    # Add a new 'Strength' column
    dplyr::select(-Pval)                          # Optionally remove the original 'Pval' column

  # Return the processed data frame
  return(expanded_grn_data)
}

