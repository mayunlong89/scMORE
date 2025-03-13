#' @title Convert scMORE TRS Results Format for Calculating Energy Score
#' @description
#' Transform scMORE TRS results for the getDistance() function by manually choosing a target cell type.
#' The TRS values of other cell types are averaged as the non-target cell type.
#' @param regulon2disease_results TRS results from regulon2disease() function.
#' @param targetCelltype Can be specified as either the index of the target cell type column (e.g., 1),
#'        or the name of the target cell type (e.g., "Monocytes").
#' @export
convertRegulon <- function(regulon2disease_results, targetCelltype = 1) {
  library(dplyr)
  library(tidyr)

  # Step 1: Extract relevant columns
  regulonScore <- regulon2disease_results[, c("RegulonName", "TRS", "Celltype")]

  # Step 2: Convert to wide format with cell types as columns
  regulonScore_convert <- regulonScore %>%
    pivot_wider(names_from = Celltype, values_from = TRS, values_fill = NA)

  # Step 3: Identify the target cell type
  celltype_pick <- colnames(regulonScore_convert)[-1]  # Exclude "RegulonName" column

  if (is.numeric(targetCelltype)) {
    if (targetCelltype < 1 || targetCelltype > length(celltype_pick)) {
      stop("The specified target cell type index is out of range.")
    }
    target_celltype <- celltype_pick[targetCelltype]
  } else if (is.character(targetCelltype)) {
    if (!(targetCelltype %in% celltype_pick)) {
      stop("The specified target cell type does not exist in the data.")
    }
    target_celltype <- targetCelltype
  } else {
    stop("Invalid targetCelltype. Must be either a numeric index or a character string.")
  }

  # Step 4: Compute Target and NonTarget scores
  non_target_columns <- setdiff(colnames(regulonScore_convert), c("RegulonName", target_celltype))

  result <- regulonScore_convert %>%
    mutate(
      Target = .[[target_celltype]],  # Extract target cell type values
      NonTarget = rowMeans(as.matrix(dplyr::select(., all_of(non_target_columns))), na.rm = TRUE)  # Compute mean for non-target cell types
    ) %>%
    dplyr::select(RegulonName, Target, NonTarget)  # Ensure correct column selection

  return(result)
}

