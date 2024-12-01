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

  # Bin sampling
  bins <- cut(
    importance_weights,
    breaks = quantile(importance_weights, probs = seq(0, 1, 0.1), na.rm = TRUE),
    include.lowest = TRUE
  )

  # Check if bins are valid
  if (any(is.na(bins))) {
    stop("Invalid bins detected. Check the importance_weights data.")
  }

  # Calculate bin probabilities
  bin_probabilities <- prop.table(table(bins))

  # Random sampling based on bin probabilities
  random_bins <- sample(
    names(bin_probabilities),
    size = length(background_genes),
    replace = TRUE,
    prob = bin_probabilities
  )

  # Generate random weights
  random_importance_weights <- sapply(random_bins, function(bin) {
    # Extract the range (lower and upper bounds) of the bin
    range <- as.numeric(unlist(regmatches(bin, gregexpr("-?[0-9.]+(?:e-?[0-9]+)?", bin, perl = TRUE))))

    # Ensure the range is valid and not too small
    if (length(range) != 2 || any(is.na(range)) || diff(range) < 1e-6) {
      stop(paste("Invalid or too small range in bin:", bin))
    }

    # Generate a random value uniformly within the bin range
    runif(1, min = range[1], max = range[2])
  })

  # Add random importance weights to random_scores
  random_scores$importance_weights_background <- random_importance_weights

  return(random_scores)
}

