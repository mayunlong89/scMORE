#' @title Convert scMORE TRS Results Format for Calculating Energy Score
#' @description
#' Transform scMORE TRS results for the getDistance() function by manually choosing a target cell type.
#' The TRS values of other cell types are averaged as the non-target cell type.
#' @param regulon2disease_results TRS results from regulon2disease() function.
#' @param targetCelltype Can be specified as either the index of the target cell type column (e.g., 1),
#'        or the name of the target cell type (e.g., "Monocytes").
#' @export
convertRegulon <- function(regulon2disease_results, targetCelltype = 1) {
  # Step 1: Extract relevant columns
  regulonScore <- regulon2disease_results[, c("RegulonName", "TRS", "Celltype")]

  # Step 2: Convert to wide format with cell types as columns
  regulonScore_convert <- regulonScore %>%
    pivot_wider(names_from = Celltype, values_from = RegulonScore, values_fill = NA)

  # Step 3: Identify the target cell type
  celltype_pick <- colnames(regulonScore_convert)[-1]  # Exclude "RegulonName" column

  if (is.numeric(targetCelltype)) {
    # If targetCelltype is an index, select by position
    target_celltype <- celltype_pick[targetCelltype]
  } else if (is.character(targetCelltype)) {
    # If targetCelltype is a name, validate and use it directly
    if (!(targetCelltype %in% celltype_pick)) {
      stop("The specified target cell type does not exist in the data.")
    }
    target_celltype <- targetCelltype
  } else {
    stop("Invalid targetCelltype. Must be either a numeric index or a character string.")
  }

  # Step 4: Calculate Target and NonTarget scores
  result <- regulonScore_convert %>%
    rowwise() %>%
    mutate(
      Target = get(target_celltype),  # Extract values for the target cell type
      NonTarget = mean(c_across(-c(RegulonName, all_of(target_celltype))), na.rm = TRUE)  # Calculate mean for non-target cell types
    ) %>%
    select(RegulonName, Target, NonTarget)  # Retain only the relevant columns

  # Step 5: Return the final result
  return(result)
}

