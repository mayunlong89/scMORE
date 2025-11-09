#' Generate Random eRegulon Scores
#'
#' This function generates random eRegulon scores for a given set of transcription factors (TFs),
#' target gene background scores, and real specificity/importance scores. It computes three scores:
#' SpecificityScore (M1, CTS), ImportanceWeightScore (M2,GRS), and the final RegulonScore (TRS), which includes a penalty term.
#'
#' @param tf_list A vector of transcription factors (TFs) to randomly select from.
#' @param target_scores_background A data frame containing background scores for target genes.
#'        The data frame must have columns `genes` and scores for filtering and sampling.
#' @param real_specificity A numeric vector of real specificity scores for sampling.
#' @param real_importance A numeric vector of real importance scores for sampling.
#' @param len_of_regulon An integer specifying the total number of genes (including the TF) in the eRegulon.
#'        Must be greater than 1.
#' @param theta A numeric value representing the weight for combining TF specificity with gene scores.
#'        Default is 0.5.
#' @param alpha A numeric value representing the weight for the penalty term in the final score.
#'        Default is 1.
#' @param top_n An integer specifying the number of top genes (ranked by importance) to include in TF calculations.
#'        Default is 5.
#'
#' @return A list containing the following elements:
#' \item{SpecificityScore}{The calculated specificity score (M1) of the regulon.}
#' \item{ImportanceWeightScore}{The calculated importance weight score (M2) of the regulon.}
#' \item{RegulonScore}{The final calculated regulon score, which incorporates a penalty term.}
#'
#' @export
#'
getRandomScore <- function(tf_list,
                           target_scores_background,
                           real_specificity,
                           real_importance,
                           len_of_regulon,
                           theta = 0.5,
                           alpha = 1,
                           top_n = 5) {

  # Step 1: Randomly select one TF
  selected_tf <- sample(tf_list, 1)

  # Step 2: Randomly select genes for the regulon, excluding the selected TF
  target_scores_background2 <- target_scores_background[target_scores_background$genes != selected_tf, ]
  len_of_regulon_genes <- len_of_regulon - 1
  if (len_of_regulon_genes <= 0) {
    stop("The length of the regulon must be greater than 1.")
  }
  sampled_genes <- sample(target_scores_background2$genes, len_of_regulon_genes, replace = TRUE)

  # Step 3: Sample random specificity and importance scores for the selected genes
  sampled_specificity <- sample(real_specificity, len_of_regulon_genes, replace = TRUE)
  sampled_importance <- sample(real_importance, len_of_regulon_genes, replace = TRUE)

  # Combine sampled TF and gene data
  sampled_data <- data.frame(
    genes = c(selected_tf, sampled_genes),
    scores = c(sample(real_specificity, 1), sampled_specificity),
    Importance_weighted = c(sample(real_importance, 1), sampled_importance),
    anno = c("TF", rep("Gene", len_of_regulon_genes))
  )

  # Step 4: Calculate scores for the TF using Top N Importance Weighted
  top_targets <- sampled_data %>%
    filter(anno == "Gene") %>%
    arrange(desc(Importance_weighted)) %>%
    slice_head(n = min(top_n, nrow(sampled_data) - 1))

  if (nrow(top_targets) < top_n) {
    warning("The number of top targets is less than the specified top_n.")
  }

  tf_specificity <- sampled_data %>% filter(anno == "TF") %>% pull(scores)
  tf_importance <- sum(top_targets$Importance_weighted, na.rm = TRUE) / nrow(top_targets)

  # Step 5: Calculate scores for the genes
  gene_scores <- sampled_data %>%
    filter(anno == "Gene") %>%
    summarize(
      m1_genes = sum(scores, na.rm = TRUE) / sqrt(n()),               # Specificity score for genes
      m2_genes = sum(Importance_weighted, na.rm = TRUE) / sqrt(n())  # Genetic risk score for genes
    )

  # Step 6: Compute M1 and M2
  M1 <- tf_specificity + theta * gene_scores$m1_genes
  M2 <- tf_importance + theta * gene_scores$m2_genes

  # Step 7: Compute the penalty term (standard deviation of M1 and M2)
  mean_gs <- (M1 + M2) / 2
  sd_mg_ms <- sqrt((M1 - mean_gs)^2 + (M2 - mean_gs)^2)

  # Step 8: Calculate the final TARS (Regulon Score)
  regulon_score <- M1 + M2 - alpha * sd_mg_ms

  # Return the results as a list
  return(list(
    SpecificityScore = M1,
    ImportanceWeightScore = M2,
    RegulonScore = regulon_score
  ))
}
