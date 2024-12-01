#'
#' @title Generate random weights for permutation
#'
#' @param regulons regulon matrix from 'snp2peak_map'
#' @param random_scores genes random Zscore and specificity score
#' @return The function gives random importance weights for background genes
#' @import fitdistrplus
#' @export
#'
getRandomWeight <- function(regulons, random_scores) {

  # Extract importance weights and background genes
  importance_weights <- regulons$Importance_weighted
  background_genes <- random_scores$genes

  # filter NA and 0 value
  importance_weights <- importance_weights[!is.na(importance_weights) & importance_weights > 0]

  # Fit a gamma distribution to importance weights
  fit_gamma <- fitdistrplus::fitdist(importance_weights, "gamma")

  # Generate random importance weights directly from the gamma distribution
  random_importance_weights <- rgamma(
    length(background_genes),
    shape = fit_gamma$estimate["shape"],
    rate = fit_gamma$estimate["rate"]
  )

  # Add random importance weights to random_scores
  random_scores$importance_weights_background <- random_importance_weights

  return(random_scores)
}
