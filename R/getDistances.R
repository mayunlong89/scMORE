#'
#' @title calcaulte Euclidean distance
#' @description
#' Define a function to calculate pairwise distances
#'
#' @param regulon_set1 regulons' scMORE score in target cell type
#' @param regulon_set2 regulons' scMORE score in other cell types
#' @export

getDistances <- function(regulon_set1, regulon_set2) {
  # Combine the two datasets (regulon_set1 and regulon_set2) and compute the distance matrix
  # `dist()` calculates pairwise distances, and `as.matrix()` converts it into a matrix
  dist_matrix <- as.matrix(dist(rbind(regulon_set1, regulon_set2)))

  # Get the number of rows (observations) in each dataset
  n1 <- nrow(regulon_set1)  # Number of points in the first dataset
  n2 <- nrow(regulon_set2)  # Number of points in the second dataset

  # Extract pairwise distances between points from the two datasets
  # These are distances between every point in `regulon_set1` and every point in `regulon_set2`
  dist_between <- dist_matrix[1:n1, (n1 + 1):(n1 + n2)]

  # Extract pairwise distances within the first dataset
  # These are distances between every pair of points in `regulon_set1`
  dist_within1 <- dist_matrix[1:n1, 1:n1]

  # Extract pairwise distances within the second dataset
  # These are distances between every pair of points in `regulon_set2`
  dist_within2 <- dist_matrix[(n1 + 1):(n1 + n2), (n1 + 1):(n1 + n2)]

  # Return a list containing all three distance matrices:
  # - `between`: Pairwise distances between `regulon_set1` and `regulon_set2`
  # - `within1`: Pairwise distances within `regulon_set1`
  # - `within2`: Pairwise distances within `regulon_set2`
  list(
    between = dist_between,
    within1 = dist_within1,
    within2 = dist_within2
  )
}
