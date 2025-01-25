#' @title Calculates the Trait-associated Regulon Score (TRS) for a given regulon
#'
#' @description
#' This function combines the importance weight (genetic risk score) and specificity score
#' (cell type specificity) to calculate a TRS (Trait-associated Regulon Score) for a given regulon.
#'
#' @param each_module_score A single regulon matrix containing importance weights and specificity scores;
#'                          this is generated in the 'regulon2disease()' function.
#' @param alternative_regulon_score_sub the cell type-specificity of each regulon calculated from the alternative methods.
#' @param theta Weight for integrating TF and gene scores (range: 0.1 to 1, default = 0.5).
#' @param alpha Flexibility in penalization (default = 1), used to scale the penalty term.
#'
#' @return A list containing:
#'         - `SpecificityScore`: The computed M1 value (specificity score).
#'         - `GeneRiskScore`: The computed M2 value (genetic risk score).
#'         - `RegulonScore`: The final TRS score for the regulon.
#'
#' @export
#'
alternativeRegulonScore <- function(each_module_score,
                                    alternative_regulon_score_sub,
                                    theta = 0.5, alpha = 1) {
  # Ensure the input data contains the required columns
  # Step 1: Extract scores for the TF
  # m1(TF) is the mean specificity score for the TF
  # m2(TF) is the mean genetic risk score for the TF

  tf_score <- each_module_score %>%
    filter(anno == "TF") %>%
    summarize(
      m2_tf = mean(Importance_weighted, na.rm = TRUE)    # m2(TF): Genetic risk score
    )

  # Step 2: Calculate scores for the genes in the regulon
  # m1(genes) and m2(genes) are calculated as the sum of scores divided by the square root of the number of genes
  gene_scores <- each_module_score %>%
    filter(anno == "Gene") %>%
    summarize(
      m2_genes = sum(Importance_weighted, na.rm = TRUE) / sqrt(n())  # Genetic risk score for genes
    )

  # Step 3: Compute M1 and M2
  # M1 combines the TF specificity score and the weighted gene specificity score
  # M2 combines the TF genetic risk score and the weighted gene genetic risk score
  M1 <- alternative_regulon_score_sub
  M2 <- tf_score$m2_tf + theta * gene_scores$m2_genes

  # Step 4: Compute the mean of M1 and M2
  mean_gs <- (M1 + M2) / 2

  # Step 5: Compute the penalty term (standard deviation of M1 and M2)
  # The penalty term adjusts for deviation between M1 and M2
  sd_mg_ms <- sqrt((M1 - mean_gs)^2 + (M2 - mean_gs)^2)

  # Step 6: Calculate the final TARS (Regulon Score)
  regulon_score <- M1 + M2 - alpha * sd_mg_ms  # Incorporate penalty scaled by alpha

  # Return the results as a list
  return(list(
    SpecificityScore = M1,        # Final M1 (specificity score)
    GeneRiskScore = M2,          # Final M2 (genetic risk score)
    RegulonScore = regulon_score # Final TRS score
  ))
}
