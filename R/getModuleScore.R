#'
#' @title Calculates the peak importance score for a given regulon
#'
#' @description
#' It extracts the highest importance score for each target and calculates the total importance for the transcription factor (TF).
#'
#' @param eachModule A single regulon; this is generated in the 'regulon2disease()' function.
#' @return A data frame with two columns: `Target` and `Importance_weighted`.
#'         The first row contains the TF name and its total importance score, followed by the targets with their highest importance scores.
#' @export
#'
getModuleScore <- function(eachModule) {
  # Ensure the input data frame contains the required columns
  if (!all(c("Target", "Importance_weighted", "TF") %in% colnames(eachModule))) {
    stop("Input data frame must contain columns: Target, Importance_weighted, and TF.")
  }

  # Extract the highest Importance_weighted value for each Target
  filtered_eachModule <- eachModule %>%
    dplyr::group_by(Target) %>%
    dplyr::slice_max(Importance_weighted, with_ties = FALSE) %>% # Ensure only the maximum value per Target is kept
    ungroup()

  # Calculate the total importance score for the TF
  # Avoid errors caused by NA values
  tf_importance <- sum(filtered_eachModule$Importance_weighted, na.rm = TRUE)/length(filtered_eachModule$Importance_weighted)

  # Extract the Target and Importance_weighted columns
  target_importance <- filtered_eachModule %>%
    dplyr::select(Target, Importance_weighted)

  # Add the TF information as the first row
  eachModule_Importance_score <- target_importance %>%
    dplyr::add_row(Target = unique(filtered_eachModule$TF), Importance_weighted = tf_importance, .before = 1)

  # Return the result
  return(eachModule_Importance_score)
}
