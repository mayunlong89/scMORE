#'
#' @title Monte Carlo permutation analysis
#'
#' @description
#' Monte Carlo permutation for calculating random specificity * JSI scores
#'
#' @param target_scores_background Data frame containing specificity scores of genes for each cell type
#' @param tf_list_1 List of high-quality TFs excluding the target TF
#' @param len_of_regulon Number of genes and TFs in the regulon
#' @param all_genes Vector of all genes from phenotype associations (e.g., MAGMA results)
#' @param top_genes Number of top-ranked disease-specific genes for JSI enrichment (default = 500)
#' @param theta Weight to adjust the specificity score of target genes in each regulon
#'              Range: 0.1 ~ 1 (default = 0.5)
#' @param pow Soft-power factor to emphasize strong specificity (default = 1)
#' @param mo Statistical model selection:
#'           - mo = 1: Genetic weight model
#'           - mo = 0: Non-genetic weight model
#'
#' @return Randomly computed specificity * JSI score for a regulon-disease link
#' @export
#'
permutation <- function(target_scores_background,
                        tf_list_1,
                        len_of_regulon,
                        all_genes,
                        top_genes = 500,
                        theta = 0.5,
                        pow = 1,
                        mo = 1) {

  # Step 1: Randomly select one TF
  selected_tf <- sample(tf_list_1, 1)

  # Step 2: Randomly select genes for the regulon, excluding the selected TF
  target_scores_background2 <- target_scores_background[target_scores_background$genes != selected_tf, ]
  len_of_regulon_genes <- len_of_regulon - 1
  sampled_genes <- sample(target_scores_background2$genes, len_of_regulon_genes)
  regulon_genes_sample <- c(selected_tf, sampled_genes)

  # Step 3: Calculate average specificity score for the sampled regulon genes
  sampled_data <- target_scores_background2[target_scores_background2$genes %in% sampled_genes, ]

  # Debugging: Check sampled data components
  #print(sampled_data$scores)
  #print(sampled_data$magma_zscore)
  #print(sampled_data$importance_weights_background)

  # Gene score calculation
  gene_scores <- (sampled_data$random_scores^pow) * (sampled_data$random_zscores * sampled_data$importance_weights_background)^mo
  #print(gene_scores)  # Debugging: Check calculated gene scores

  avg_gene_score <- mean(gene_scores, na.rm = TRUE)

  # Step 4: Calculate specificity score for the selected TF
  tf_data <- target_scores_background[target_scores_background$genes == selected_tf, c("random_scores", "random_zscores")]
  tf_importance_weights_background <- sum(sampled_data$importance_weights_background, na.rm = TRUE)/length(sampled_data$importance_weights_background)
  tf_score <- as.numeric((tf_data$random_scores^pow) * (tf_data$random_zscores* tf_importance_weights_background)^mo)
  #tf_score <- as.numeric((tf_data$scores^pow) * (tf_data$magma_zscore)^mo)

  # Step 5: Compute overall regulon specificity score
  regulon_score <- tf_score + theta * avg_gene_score

  # Step 6: Randomly select background genes for JSI calculation
  bg_genes <- sample(all_genes, top_genes)
  inter_genes <- length(intersect(bg_genes, regulon_genes_sample))
  union_genes <- length(union(bg_genes, regulon_genes_sample))
  JSI_sample <- (inter_genes + 1) / (union_genes + 1)

  # Step 7: Compute specificity * JSI score for the regulon-disease link
  ct_score_sample <- regulon_score * JSI_sample
  #ct_score_sample <- regulon_score
  # Return the computed score
  return(ct_score_sample)
}

