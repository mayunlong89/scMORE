#'
#' @title Identifying Cell Type-Specific Regulons Relevant to Disease
#'
#' @description
#' This function identifies cell type-specific TF-regulons relevant to a disease using MAGMA-based gene risk scores, peak-gene and SNP-peak mapping, and COSR and JSI scoring.
#'
#' @param grn_outputs Output from GRN analysis containing TF names and regulatory connections.
#' @param target_scores Specificity scores for genes and TFs across cell types.
#' @param snp_info SNP information used for linking SNPs to peaks.
#' @param geneRiskScores MAGMA-based gene risk scores with columns SYMBOL, logP, and ZSTAT.
#' @param perm_n Number of permutations for Monte Carlo analysis (default = 1000).
#' @param top_genes Number of top-ranked genes to consider for regulon-disease links (default = 500).
#' @param theta Weighting factor for averaging scores (range: 0.1 to 1, default = 0.5).
#' @param pow Power factor to amplify specificity differences (default = 1).
#' @param mo Statistical model: 1 = genetic weight model, 0 = non-genetic weight model (default = 1).
#' @param buffer Buffer size for SNP-peak mapping (default = 50 bp).
#' @return A list containing regulon-disease link scores (MORE_score) and their Monte Carlo p-values (MORE_perm_Pval).
#' @export
#'
regulon2disease <- function(grn_outputs, target_scores, snp_info, geneRiskScores,
                            perm_n = 1000, top_genes = 500, theta = 0.5, pow = 1, mo = 1, buffer = 50) {

  # Step 1: Map SNPs to TF-peaks and peaks to genes
  peak2gene_strength <- peak2gene(grn_outputs)  # Map peaks to genes based on GRN outputs
  snp2peak_map <- snp2peak(snp_info, peak2gene_strength, buffer = buffer)  # Map SNPs to peaks
  snp2peak_map <- getPeakScore(snp2peak_map)  # Calculate peak importance scores

  # Step 2: Extract regulons
  regulons <- snp2peak_map[, c("snp_id", "peak_ids", "TF", "Target", "Importance_weighted")]
  snp2peak_map <- NULL  # Free memory
  peak2gene_strength <- NULL

  # Step 3: Initialize results data frames
  tf_list <- grn_outputs$tf_names
  regulon_names <- data.frame(paste(tf_list, "_regulon", sep = ""))
  Final_regulon_score <- data.frame(ID_regulon = tf_list)
  Final_regulon_MORE_score <- data.frame(ID_regulon = tf_list)
  Final_regulon_MORE_perm_p <- data.frame(ID_regulon = tf_list)

  # Step 4: Normalize target scores across cell types
  all_celltype_names <- unique(target_scores[, 3])  # Get all cell type names
  for (m in all_celltype_names) {
    target_scores$scores[target_scores$celltypes == m] <- max_min_scale(target_scores$scores[target_scores$celltypes == m])
  }

  # Step 5: Initialize progress bar
  pb <- txtProgressBar(style = 3)
  total_run <- length(all_celltype_names) * length(tf_list)
  count <- 0
  start_time <- Sys.time()

  # Step 6: Calculate regulon specificity and disease relevance
  for (cell_type in all_celltype_names) {
    regulon_MORE_score <- numeric()
    regulon_MORE_perm_p <- numeric()
    regulon_score_all <- numeric()

    for (tf in tf_list) {
      # Extract genes in the regulon
      Module <- regulons[regulons$TF == tf, ]
      Module_regulon <- unique(c(Module$TF, Module$Target))

      # Calculate importance scores for the regulon
      eachModule_Importance_score <- getModuleScore(Module)

      # Filter and annotate target scores with MAGMA z-scores
      target_scores_sub <- target_scores[target_scores$celltypes == cell_type, ]
      geneRiskScores_sub <- geneRiskScores[geneRiskScores$SYMBOL %in% target_scores_sub$genes, ]
      target_scores_sub <- target_scores_sub[target_scores_sub$genes %in% geneRiskScores_sub$SYMBOL, ]
      target_scores_sub$magma_zscore <- geneRiskScores_sub$ZSTAT[match(target_scores_sub$genes, geneRiskScores_sub$SYMBOL)]

      # Combine scores for the TF and target genes
      each_module_score <- target_scores_sub[match(target_scores_sub$genes, Module_regulon), ]
      each_module_score$Importance_weighted <- eachModule_Importance_score$Importance_weighted[match(each_module_score$genes, eachModule_Importance_score$Target)]

      # Calculate combined scores and regulon specificity
      tf_combined_score <- calculate_tf_combined_score(each_module_score, theta, pow, mo)
      regulon_score <- calculate_regulon_score(tf_combined_score, each_module_score, theta, pow, mo)

      # Calculate JSI and MORE score
      top_ranked_genes <- geneRiskScores$SYMBOL[1:top_genes]
      JSI_score <- calculate_JSI(top_ranked_genes, Module_regulon)
      MORE_score <- regulon_score * JSI_score

      # Monte Carlo permutation
      perm_results <- replicate(perm_n, permutation(target_scores, tf_list, Module_regulon, all_genes, top_genes, theta, pow, mo))
      perm_p <- calculate_perm_p(MORE_score, perm_results)

      # Save scores
      regulon_MORE_score <- c(regulon_MORE_score, MORE_score)
      regulon_MORE_perm_p <- c(regulon_MORE_perm_p, perm_p)

      # Update progress bar
      count <- count + 1
      setTxtProgressBar(pb, count / total_run)
    }

    # Save results for this cell type
    save_results(cell_type, regulon_MORE_score, regulon_MORE_perm_p, Final_regulon_score, Final_regulon_MORE_score, Final_regulon_MORE_perm_p)
  }

  # Finalize and return results
  close(pb)
  end_time <- Sys.time()
  print(paste("Running time:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds"))

  list(MORE_score = Final_regulon_MORE_score, MORE_perm_Pval = Final_regulon_MORE_perm_p)
}
