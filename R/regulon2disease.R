#' @title Identify cell type-specific regulons relevant to disease
#'
#' @param grn_outputs GRN outputs containing TF-gene relationships
#' @param target_scores Matrix of genes and TF specificity scores across cell types (columns: genes, scores, celltypes)
#' @param snp_info Information about SNPs for the analysis
#' @param geneRiskScores MAGMA-based gene association results (columns: SYMBOL, logP, ZSTAT)
#' @param perm_n Number of permutations for Monte Carlo simulation (default = 1000)
#' @param theta Weight for integrating TF and gene scores (range: 0.1~1, default = 0.5)
#' @param alpha Flexibility in penalization (default = 1)
#' @param top_n Top n targets of each TF to calculate the importance of the TF (default = 5)
#' @param buffer Distance buffer (in base pairs) for SNP-to-peak mapping (default = 500bp), which means each peak will be extended by 500bp upstream and 500bp downstream.
#' @return A data frame containing specificity, genetic risk score and regulon score and corresponding P values.
#' @export
#'
regulon2disease <- function(grn_outputs,
                            target_scores,
                            geneRiskScores,
                            snp_info,
                            perm_n = 1000,
                            theta = 0.5,
                            alpha = 1,
                            top_n = 5,
                            buffer = 500) {

  # Step 1: Map SNPs to TF-peaks and target genes
  message("Step 1: Mapping SNPs to TF-peaks and target genes...")
  peak2gene_strength <- peak2gene(grn_outputs)  # Map peaks to genes
  snp2peak_map <- snp2peak(snp_info, peak2gene_strength, buffer = buffer)  # SNP-to-peak mapping

  # Add gene risk scores to SNP-to-peak mapping
  snp2peak_map$geneScores <- geneRiskScores$logP[match(snp2peak_map$Target, geneRiskScores$SYMBOL)]
  snp2peak_map <- getPeakScore(snp2peak_map)  # Calculate peak importance scores

  # Extract regulons containing SNP, peak, TF, target gene, and importance score
  regulons <- snp2peak_map[, c("snp_id", "peak_ids", "TF", "Target", "Importance_weighted")]

  # Step 2: Define TF list and cell types
  message("Step 2: Defining TF list and cell types...")
  tf_list <- grn_outputs$tf_names  # Transcription factors
  all_celltype_names <- unique(target_scores[, "celltypes"])  # Cell types

  # Initialize results data frame
  all_regulon_results_df <- data.frame(
    RegulonID = character(),
    RegulonName = character(),
    SpecificityScore = numeric(),
    SpecificityScore_p = numeric(),
    GeneRiskScore = numeric(),
    ImportanceWeightScore_p = numeric(),
    RegulonScore = numeric(),
    RegulonScore_p = numeric(),
    Celltype = character(),
    stringsAsFactors = FALSE
  )

  # Open progress bar
  pb <- txtProgressBar(style = 3)
  start_time <- Sys.time()
  total_run <- length(all_celltype_names) * length(tf_list)
  count <- 0

  # Step 3: Perform COSR and MORE analysis for each cell type and TF
  for (i in seq_along(all_celltype_names)) {
    for (j in seq_along(tf_list)) {

      # Extract regulon genes for the current TF
      Module <- regulons[regulons$TF == tf_list[j], ]
      Module_regulon <- c(unique(Module$TF), unique(Module$Target))

      # Calculate importance scores for the regulon
      eachModule_Importance_score <- getRiskScore(Module, top_n = top_n)

      # Subset target scores for the current cell type
      target_scores_sub <- target_scores[target_scores[, "celltypes"] == all_celltype_names[i], ]

      # Filter specificity scores for regulon genes
      each_module_score <- target_scores_sub[match(target_scores_sub[, "genes"], Module_regulon, nomatch = 0) > 0, ]
      each_module_score$anno <- ifelse(each_module_score$genes == Module_regulon[1], "TF", "Gene")

      # Add importance scores to regulon
      each_module_score$Importance_weighted <- eachModule_Importance_score$Importance_weighted[match(
        each_module_score$genes,
        eachModule_Importance_score$Target
      )]

      # Compute regulon scores (TARS)
      TARS <- getRegulonScore(each_module_score, theta = theta, alpha = alpha)

      # Run Monte Carlo permutation analysis
      len_of_regulon <- length(Module_regulon)
      target_scores_background <- target_scores_sub
      real_specificity <- target_scores_sub$scores[match(target_scores_sub[, "genes"], regulons$Target, nomatch = 0) > 0]
      real_importance <- regulons$Importance_weighted

      perm_results <- replicate(
        perm_n,
        getRandomScore(
          tf_list,
          target_scores_background,
          real_specificity,
          real_importance,
          len_of_regulon,
          theta = theta,
          alpha = alpha,
          top_n = top_n
        )
      )

      # Extract and clean permutation scores
      perm_specificity <- as.numeric(perm_results["SpecificityScore", ])
      perm_importance <- as.numeric(perm_results["ImportanceWeightScore", ])
      perm_regulon <- as.numeric(perm_results["RegulonScore", ])

      # Remove NA values
      perm_specificity <- perm_specificity[!is.na(perm_specificity)]
      perm_importance <- perm_importance[!is.na(perm_importance)]
      perm_regulon <- perm_regulon[!is.na(perm_regulon)]

      # Calculate p-values and z-scores
      p_specificity <- (1 + sum(perm_specificity >= TARS$SpecificityScore)) / (1 + length(perm_specificity))
      p_importance <- (1 + sum(perm_importance >= TARS$GeneRiskScore)) / (1 + length(perm_importance))
      p_regulon <- (1 + sum(perm_regulon >= TARS$RegulonScore)) / (1 + length(perm_regulon))

      z_specificity <- if (sd(perm_specificity) == 0) NA else (TARS$SpecificityScore - mean(perm_specificity)) / sd(perm_specificity)
      z_importance <- if (sd(perm_importance) == 0) NA else (TARS$GeneRiskScore - mean(perm_importance)) / sd(perm_importance)
      z_regulon <- if (sd(perm_regulon) == 0) NA else (TARS$RegulonScore - mean(perm_regulon)) / sd(perm_regulon)

      # Append results to the data frame
      all_regulon_results_df <- rbind(
        all_regulon_results_df,
        data.frame(
          RegulonID = paste0("Regulon_", j),
          RegulonName = Module_regulon[1],
          SpecificityScore = z_specificity,
          SpecificityScore_p = p_specificity,
          GeneRiskScore = z_importance,
          ImportanceWeightScore_p = p_importance,
          RegulonScore = z_regulon,
          RegulonScore_p = p_regulon,
          Celltype = all_celltype_names[i]
        )
      )

      # Update progress bar
      count <- count + 1
      setTxtProgressBar(pb, count / total_run)
    }
  }

  # Close progress bar and record running time
  close(pb)
  end_time <- Sys.time()
  cat("Running time:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")

  # Add significance column
  all_regulon_results_df$Significance <- ifelse(
    all_regulon_results_df$SpecificityScore_p < 0.05 &
      all_regulon_results_df$ImportanceWeightScore_p < 0.05 &
      all_regulon_results_df$RegulonScore_p < 0.05,
    "Significant",
    "Nonsignificant"
  )

  return(all_regulon_results_df)
}
