#'
#' @title Calculate the random scores for background genes
#'
#' @description
#' Generate random scores based on the fitted distributions of specificity scores and MAGMA z-scores.
#'
#' @param target_scores_sub A data frame containing specificity scores and MAGMA z-scores for each gene.
#'                          Must contain the following columns:
#'                          - scores: Specificity scores for each gene.
#'                          - magma_zscore: MAGMA z-scores for each gene.
#' @return A data frame with the original gene names and generated random specificity and MAGMA z-scores.
#' @export
#'
getRandomScore <- function(target_scores_sub) {
  # Step 1: Fit the distribution of specificity scores (Gamma distribution)
  fit_scores <- fitdist(target_scores_sub$scores, "gamma")
  #summary(fit_scores)

  # Step 2: Fit the distribution of MAGMA z-scores (Normal distribution)
  fit_zscore <- fitdist(target_scores_sub$magma_zscore, "norm")
  #summary(fit_zscore)

  # Step 3: Generate random data based on the fitted distributions
  set.seed(1234)  # Set random seed for reproducibility

  random_scores <- rgamma(
    n = nrow(target_scores_sub),
    shape = fit_scores$estimate["shape"],
    rate = fit_scores$estimate["rate"]
  )

  random_zscores <- rnorm(
    n = nrow(target_scores_sub),
    mean = fit_zscore$estimate["mean"],
    sd = fit_zscore$estimate["sd"]
  )

  # Step 4: Combine the generated random data into a data frame
  random_data <- data.frame(
    genes = target_scores_sub$genes,
    random_scores = random_scores,
    random_zscores = random_zscores
  )

  # Return the result
  return(random_data)
}
