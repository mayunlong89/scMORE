
#'
#' @title Energy score calculation
#'
#' @description
#' Calculate the energy statistics of regulons in target cell type compared to other cell types
#'
#' @param regulon_set1 regulons' scMORE score in target cell type
#' @param regulon_set2 regulons' scMORE score in other cell types
#' @export

getEnergyScore <- function(regulon_set1, regulon_set2) {
  # Calculate pairwise distances between and within the datasets regulon_set1 and regulon_set2
  distances <- getDistances(regulon_set1, regulon_set2)

  # Compute δ_XY: Mean of the pairwise distances between points in X and points in Y
  delta_set12 <- mean(distances$between)

  # Compute σ_X: Mean of the pairwise distances within dataset X (bias-corrected)
  sigma_set1 <- sum(distances$within1[lower.tri(distances$within1)]) / (nrow(regulon_set1) * (nrow(regulon_set1) - 1))

  # Compute σ_Y: Mean of the pairwise distances within dataset Y (bias-corrected)
  sigma_set2 <- sum(distances$within2[lower.tri(distances$within2)]) / (nrow(regulon_set2) * (nrow(regulon_set2) - 1))

  # Apply the energy distance formula: E(X, Y) = 2 * δ_XY - σ_X - σ_Y
  Escore <- 2 * delta_set12 - sigma_set1 - sigma_set2

  return(Escore)
}
