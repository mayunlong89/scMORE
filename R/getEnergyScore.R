#' @title Energy Score Calculation
#' @description
#' Calculate the energy statistics of regulons in the target cell type compared to other cell types.
#' @param regulon2disease_results TRS results from regulon2disease() function.
#' @param targetCelltype Can be specified as either the index of the target cell type column (e.g., 1),
#'        or the name of the target cell type (e.g., "Monocytes").
#' @export
#'
getEnergyScore <- function(regulon2disease_results, targetCelltype = 1) {

  # Convert format
  target_results <- convertRegulon(regulon2disease_results, targetCelltype = targetCelltype)

  # Calculate pairwise distances between and within the datasets regulon_set1 and regulon_set2
  distances <- getDistances(target_results)

  # Compute δ_XY: Mean of the pairwise distances between points in X and points in Y
  delta_set12 <- mean(distances$between, na.rm = TRUE)

  # Compute σ_X: Mean of the pairwise distances within dataset X (bias-corrected)
  n1 <- nrow(target_results)
  sigma_set1 <- sum(distances$within_target[lower.tri(distances$within_target)], na.rm = TRUE) /
    (n1 * (n1 - 1) / 2)

  # Compute σ_Y: Mean of the pairwise distances within dataset Y (bias-corrected)
  n2 <- nrow(target_results)
  sigma_set2 <- sum(distances$within_nontarget[lower.tri(distances$within_nontarget)], na.rm = TRUE) /
    (n2 * (n2 - 1) / 2)

  # Apply the energy distance formula: E(X, Y) = 2 * δ_XY - σ_X - σ_Y
  Escore <- 2 * delta_set12 - sigma_set1 - sigma_set2

  return(Escore)
}
