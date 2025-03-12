#' @title Identify Cell Type-Specific Regulons Relevant to Disease
#'
#' @description
#' Function to identify cell type-specific regulons (TRS) that are associated with disease by integrating
#' cell type specificity scores(CTS) and genetic relevance scores (GRS).
#'
#' @param grn_outputs GRN (Gene Regulatory Network) outputs containing TF-gene relationships.
#' @param target_scores A matrix of genes and their TF specificity scores across cell types.
#'                      Columns include:
#'                      - `genes`: Gene names
#'                      - `scores`: TF specificity scores
#'                      - `celltypes`: Cell type annotations
#' @param snp_info SNP information used for the analysis.
#' @param geneRiskScores Gene-based genetic association results from MAGMA. Must include the following columns:
#'                       - `SYMBOL`: Gene symbols
#'                       - `logP`: Log-transformed p-values
#'                       - `ZSTAT`: Z-scores from MAGMA
#' @param perm_n Number of permutations for Monte Carlo simulation. Default: 1000.
#' @param theta Weighting factor for integrating TF and gene specificity scores. Range: 0.1 to 1. Default: 0.5.
#' @param alpha Flexibility parameter for penalization in the scoring model. Default: 1.
#' @param top_n Number of top targets for each TF used to calculate TF importance. Default: 5.
#' @param buffer Numeric value specifying the flanking region size for genomic peaks.
#'               Extends the peak range upstream and downstream by the specified value.
#'               For example, `buffer = 500` adds 500 bp on both sides of a peak. Default: 500 bp.
#' @param infer_method GRN inference method:
#'   - 'glm' (default): Generalized linear model
#'   - 'cv.glmnet': Regularized GLMs
#'   - 'brms': Bayesian regression models
#'   - 'bagging_ridge': Bagging ridge and Bayesian ridge (CellOracle)
#'   - 'xgb': XGBoost gradient-boosted Random Forest (SCENIC)
#' @param p1 Threshold for statistical significance of the cell type-specificity score (CTS)
#'           for each regulon. Default: 0.05.
#' @param p2 Threshold for statistical significance of the genetic relevance score (GRS)
#'           for each regulon. Default: 0.05.
#' @param p3 Threshold for statistical significance of the trait-associated regulon score (TRS).
#'           Default: 0.05.
#'
#' @return A data frame containing the following:
#'         - `specificity`: Cell type-specificity scores for each regulon (CTS).
#'         - `genetic_risk_score`: Genetic relevance scores for each regulon (GRS).
#'         - `regulon_score`: Trait-associated regulon scores (TRS).
#'         - `p_values`: Statistical significance values for specificity, risk, and regulon scores.
#'
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
                            buffer = 500,
                            infer_method = 'glm',
                            p1 = 0.05,
                            p2 = 0.05,
                            p3 = 0.05) {

  # Step 4.1: Map SNPs to TF-peaks and target genes
  message("Step 4.1: Mapping SNPs to TF-peaks and target genes...")
  peak2gene_strength <- peak2gene(grn_outputs,infer_method = infer_method)  # Map peaks to genes
  snp2peak_map <- snp2peak(snp_info, peak2gene_strength, buffer = buffer)  # SNP-to-peak mapping

  # Add gene risk scores to SNP-to-peak mapping
  snp2peak_map$geneScores <- geneRiskScores$logP[match(snp2peak_map$Target, geneRiskScores$SYMBOL)]
  snp2peak_map <- getPeakScore(snp2peak_map)  # Calculate peak importance scores

  # Extract regulons containing SNP, peak, TF, target gene, and importance score (GRS score of each node)
  regulons <- snp2peak_map[, c("snp_id", "peak_ids", "TF", "Target", "Importance_weighted")]

  # Step 4.2: Define Regulon list and cell types
  message("Step 4.2: Defining Regulon list and cell types...")
  tf_list <- grn_outputs$tf_names  # Transcription factors
  all_celltype_names <- unique(target_scores[, "celltypes"])  # Cell types

  # Initialize results data frame
  all_regulon_results_df <- data.frame(
    RegulonID = character(),
    RegulonName = character(),
    CTS = numeric(),
    CTS_P_value = numeric(),
    GRS = numeric(),
    GRS_P_value = numeric(),
    TRS = numeric(),
    TRS_P_value = numeric(),
    Celltype = character(),
    stringsAsFactors = FALSE
  )

  # Step 4.3: Perform scMORE analysis for each cell type and regulon
  message("Step 4.3: Calculate Trait-Asssociated Regulon Score and Signifiance...")

  # Open progress bar
  pb <- txtProgressBar(style = 3)
  start_time <- Sys.time()
  total_run <- length(all_celltype_names) * length(tf_list)
  count <- 0

  for (i in seq_along(all_celltype_names)) {
    for (j in seq_along(tf_list)) {

      # Extract regulon genes for the current TF
      Module <- regulons[regulons$TF == tf_list[j], ]
      Module_regulon <- c(unique(Module$TF), unique(Module$Target))

      # Calculate importance scores for the regulon
      eachModule_Importance_score <- getRiskScore(Module, top_n = top_n)

      # Subset target cell type-specificity scores for the current cell type
      target_scores_sub <- target_scores[target_scores[, "celltypes"] == all_celltype_names[i], ]

      # Filter specificity scores for regulon genes
      each_module_score <- target_scores_sub[match(target_scores_sub[, "genes"], Module_regulon, nomatch = 0) > 0, ]
      each_module_score$anno <- ifelse(each_module_score$genes == Module_regulon[1], "TF", "Gene")

      # Add importance scores to regulon
      each_module_score$Importance_weighted <- eachModule_Importance_score$Importance_weighted[match(
        each_module_score$genes,
        eachModule_Importance_score$Target
      )]

      #Here, header of each_module_score: genes, scores, celltypes, anno, Importance_weighted
      #scores: cell type-specificity (CTS) score of each gene (node)
      #Importance_weighted: genetic relevance scores (GRS) of each node

      # Compute regulon scores (TRS)
      TRS <- getRegulonScore(each_module_score, theta = theta, alpha = alpha)

      # Run Monte Carlo (MC) permutation analysis
      len_of_regulon <- length(Module_regulon)
      target_scores_background <- target_scores_sub
      real_specificity <- target_scores_sub$scores[match(target_scores_sub[, "genes"], regulons$Target, nomatch = 0) > 0]
      real_importance <- regulons$Importance_weighted

      # MC analysis
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
      p_specificity <- (1 + sum(perm_specificity >= TRS$SpecificityScore)) / (1 + length(perm_specificity))
      p_importance <- (1 + sum(perm_importance >= TRS$GeneRiskScore)) / (1 + length(perm_importance))
      p_regulon <- (1 + sum(perm_regulon >= TRS$RegulonScore)) / (1 + length(perm_regulon))

      # Compute z-scores for each regulon
      z_specificity <- if (sd(perm_specificity) == 0) NA else (TRS$SpecificityScore - mean(perm_specificity)) / sd(perm_specificity)
      z_importance <- if (sd(perm_importance) == 0) NA else (TRS$GeneRiskScore - mean(perm_importance)) / sd(perm_importance)
      z_regulon <- if (sd(perm_regulon) == 0) NA else (TRS$RegulonScore - mean(perm_regulon)) / sd(perm_regulon)

      # Append results to the data frame
      all_regulon_results_df <- rbind(
        all_regulon_results_df,
        data.frame(
          RegulonID = paste0("Regulon_", j),
          RegulonName = Module_regulon[1],
          CTS = z_specificity,
          CTS_P_value = p_specificity,
          GRS = z_importance,
          GRS_P_value = p_importance,
          TRS = z_regulon,
          TRS_P_value = p_regulon,
          Celltype = all_celltype_names[i]
        )
      )

      # Update progress bar
      count <- count + 1
      setTxtProgressBar(pb, count / total_run)
    }
  }

  # Step 4.4 Add significance column
  message("Step 4.4: Adding Significance column...")
  all_regulon_results_df$Significance <- ifelse(
    all_regulon_results_df$CTS_P_value  < p1 &
      all_regulon_results_df$GRS_P_value < p2 &
      all_regulon_results_df$TRS_P_value < p3,
    "Significant",
    "Nonsignificant"
  )

  # Close progress bar and record running time
  close(pb)
  end_time <- Sys.time()
  cat("Running time:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")

  return(all_regulon_results_df)
}
