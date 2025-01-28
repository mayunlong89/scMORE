#'
#' @title Max-min scaling
#'
#' @description
#' Max-min normalization with epsilon adjustment
#'
#' @param x Numeric vector to be scaled
#' @param epsilon  Define a small epsilon to avoid zero division (default = 1e-6)
#' @return Scaled vector with values ranging from 0 to 1, adjusted to avoid all-zero results
#' @export
#'
max_min_scale <- function(x, epsilon=1e-6) {
  # Replace NA values with the mean of non-NA elements
  x[is.na(x)] <- mean(x, na.rm = TRUE)

  # Apply max-min scaling with epsilon adjustment
  scale_vec <- (x - min(x) + epsilon) / (max(x) - min(x) + epsilon)

  return(scale_vec)
}


