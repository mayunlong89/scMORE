#' Generate Random Regulon Scores
#'
#' This function generates random regulon scores for a given set of transcription factors (TFs),
#' target gene background scores, and real specificity/importance scores. It computes three scores:
#' SpecificityScore (M1, CTS), ImportanceWeightScore (M2,GRS), and the final RegulonScore (TRS), which includes a penalty term.
#'
#' @param tf_list A vector of transcription factors (TFs) to randomly select from.
#' @param Celltype the cell type used for calculating the regulon specificity.
#' @param background_genes All genes in single-cell data
#' @param real_importance A numeric vector of real importance scores for sampling.
#' @param len_of_regulon An integer specifying the total number of genes (including the TF) in the regulon.
#'        Must be greater than 1.
#' @param alternative A string specifying the method to use for scoring. Options are:
#'                "AUCell", "RSS", "VAM", "Seurat", "UCell". Default is "AUCell".
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
alternativeRandomScore <- function(tf_list,
                                   background_genes,
                                   real_importance,
                                   len_of_regulon,
                                   Celltype,
                                   alternative,
                                   theta = 0.5,
                                   alpha = 1,
                                   top_n = 5) {

  # Step 1: Randomly select one TF
  selected_tf <- sample(tf_list, 1)

  # Step 2: Randomly select genes for the regulon, excluding the selected TF
  len_of_regulon_genes <- len_of_regulon - 1
  if (len_of_regulon_genes <= 0) {
    stop("The length of the regulon must be greater than 1.")
  }
  sampled_genes <- sample(background_genes, len_of_regulon_genes, replace = TRUE)

  # Step 3: Sample random specificity and importance scores for the selected genes
  #sampled_specificity <- sample(real_specificity, len_of_regulon_genes, replace = TRUE)
  sampled_importance <- sample(real_importance, len_of_regulon_genes, replace = TRUE)

  # Combine sampled TF and gene data
  sampled_data <- data.frame(
    genes = c(selected_tf, sampled_genes),
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

  tf_importance <- sum(top_targets$Importance_weighted, na.rm = TRUE) / nrow(top_targets)

  # Step 5: Calculate scores for the genes
  gene_scores <- sampled_data %>%
    filter(anno == "Gene") %>%
    summarize(
      m2_genes = sum(Importance_weighted, na.rm = TRUE) / sqrt(n())  # Genetic risk score for genes
    )

  Module_regulon = c(selected_tf, sampled_genes)
  specificity_score <- alternativeRandomSpecificity(single_cell = single_cell,
                                                    Module_regulon = Module_regulon,
                                                    alternative = alternative)


  # Step 6: Compute M1 and M2
  M1 <-  specificity_score$scores[which(specificity_score$celltypes==Celltype)]
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
