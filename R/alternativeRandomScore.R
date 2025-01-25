#' @title Generate Random Regulon Scores using Alternative Methods
#' @description
#' This function generates random regulon scores for a given set of transcription factors (TFs),
#' background genes, and real importance scores using alternative methods.
#' It computes three scores:
#' - SpecificityScore (M1, CTS)
#' - ImportanceWeightScore (M2, GRS)
#' - RegulonScore (TRS), which includes a penalty term.
#'
#' @param tf_list A vector of transcription factors (TFs) to randomly select from.
#' @param background_genes A vector of all genes in the single-cell data.
#' @param real_importance A numeric vector of real importance scores for sampling.
#' @param len_of_regulon An integer specifying the total number of genes (including the TF) in the regulon.
#'        Must be greater than 1.
#' @param Celltype The cell type used for calculating the regulon specificity.
#' @param alternative A string specifying the method to use for scoring. Options: "AUCell", "RSS", "VAM", "Seurat", "UCell".
#'        Default: "AUCell".
#' @param theta A numeric value representing the weight for combining TF specificity with gene scores. Default: 0.5.
#' @param alpha A numeric value representing the weight for the penalty term in the final score. Default: 1.
#' @param top_n An integer specifying the number of top genes (ranked by importance) to include in TF calculations. Default: 5.
#'
#' @return A list containing:
#' - SpecificityScore: The calculated specificity score (M1) of the regulon.
#' - ImportanceWeightScore: The calculated importance weight score (M2) of the regulon.
#' - RegulonScore: The final calculated regulon score, incorporating a penalty term.
#'
#' @export
alternativeRandomScore <- function(tf_list,
                                   background_genes,
                                   real_importance,
                                   len_of_regulon,
                                   Celltype,
                                   alternative = "AUCell",
                                   theta = 0.5,
                                   j=j,
                                   alpha = 1,
                                   top_n = 5) {
  # Step 1: Randomly select one TF
  selected_tf <- sample(tf_list, 1)

  # Step 2: Randomly select genes for the regulon, excluding the TF
  len_of_regulon_genes <- len_of_regulon - 1
  if (len_of_regulon_genes <= 0) stop("The length of the regulon must be greater than 1.")
  sampled_genes <- sample(background_genes, len_of_regulon_genes, replace = TRUE)

  # Step 3: Sample random importance scores for the selected genes
  sampled_importance <- sample(real_importance, len_of_regulon_genes, replace = TRUE)

  # Combine sampled TF and gene data
  sampled_data <- data.frame(
    genes = c(selected_tf, sampled_genes),
    Importance_weighted = c(sample(real_importance, 1), sampled_importance),
    anno = c("TF", rep("Gene", len_of_regulon_genes))
  )

  # Step 4: Calculate TF scores using top N genes ranked by importance
  top_targets <- sampled_data %>%
    filter(anno == "Gene") %>%
    arrange(desc(Importance_weighted)) %>%
    slice_head(n = min(top_n, nrow(sampled_data) - 1))

  if (nrow(top_targets) < top_n) warning("The number of top targets is less than the specified top_n.")

  tf_importance <- sum(top_targets$Importance_weighted, na.rm = TRUE) / nrow(top_targets)

  # Step 5: Calculate gene scores
  gene_scores <- sampled_data %>%
    filter(anno == "Gene") %>%
    summarize(
      m2_genes = sum(Importance_weighted, na.rm = TRUE) / sqrt(n()) # Genetic risk score for genes
    )

  # Step 6: Calculate specificity score using the alternative method
  Module_regulon1 <- c(selected_tf, sampled_genes)
  specificity_score <- alternativeRandomSpecificity(
    single_cell = single_cell,
    Module_regulon = Module_regulon1,
    alternative = alternative,
    tf_list,
    j=j
  )
  #specificity_score
  # Step 7: Compute M1 (Specificity Score) and M2 (Importance Score)
  M1 <- specificity_score$scores[which(specificity_score$celltypes == Celltype)]
  M2 <- tf_importance + theta * gene_scores$m2_genes

  # Step 8: Compute the penalty term (standard deviation of M1 and M2)
  mean_gs <- (M1 + M2) / 2
  sd_mg_ms <- sqrt((M1 - mean_gs)^2 + (M2 - mean_gs)^2)

  # Step 9: Calculate the final Regulon Score
  regulon_score <- M1 + M2 - alpha * sd_mg_ms

  # Return the results as a list
  return(list(
    SpecificityScore = M1,
    ImportanceWeightScore = M2,
    RegulonScore = regulon_score
  ))
}
