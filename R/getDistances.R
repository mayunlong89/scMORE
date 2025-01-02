#' @title Compute Euclidean Distances Between Target and Non-Target Scores
#' @description
#' This function calculates the Euclidean distances between Target and NonTarget scores
#' from the input dataset, as well as the within distances for each group.
#' @param target_results A data frame with three columns: RegulonName, Target, and NonTarget.
#' @return A list containing:
#'         - `between`: Pairwise distances between Target and NonTarget.
#'         - `within_target`: Pairwise distances within the Target group.
#'         - `within_nontarget`: Pairwise distances within the NonTarget group.
#' @export
getDistances <- function(target_results) {
  # Step 1: Separate Target and NonTarget scores
  regulon_set1 <- data.frame(value = target_results$Target)
  regulon_set2 <- data.frame(value = target_results$NonTarget)

  # Step 2: Combine the two datasets for distance calculation
  combined_data <- rbind(regulon_set1, regulon_set2)

  # Step 3: Compute the distance matrix
  dist_matrix <- as.matrix(dist(combined_data))

  # Step 4: Get the number of rows in each dataset
  n1 <- nrow(regulon_set1)  # Number of points in the Target group
  n2 <- nrow(regulon_set2)  # Number of points in the NonTarget group

  # Step 5: Extract pairwise distances
  dist_between <- dist_matrix[1:n1, (n1 + 1):(n1 + n2)]  # Between Target and NonTarget
  dist_within_target <- dist_matrix[1:n1, 1:n1]          # Within Target group
  dist_within_nontarget <- dist_matrix[(n1 + 1):(n1 + n2), (n1 + 1):(n1 + n2)]  # Within NonTarget group

  # Step 6: Return the distances as a list
  list(
    between = dist_between,
    within_target = dist_within_target,
    within_nontarget = dist_within_nontarget
  )
}
