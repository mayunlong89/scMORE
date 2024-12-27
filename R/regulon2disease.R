#' @title Identify cell type-specific regulons relevant to disease using COSR and MORE methods
#'
#' @param grn_outputs GRN outputs containing TF-gene relationships
#' @param target_scores Matrix of genes and TF specificity scores across cell types (columns: genes, scores, celltypes)
#' @param snp_info Information about SNPs for the analysis
#' @param geneRiskScores MAGMA-based gene association results (columns: SYMBOL, logP, ZSTAT)
#' @param perm_n Number of permutations for Monte Carlo simulation (default = 1000)
#' @param theta Weight for integrating TF and gene scores (range: 0.1~1, default = 0.5)
#' @param alpha Flexibility in penalization (default=1) 
#' @param top_n Top n targets of each TF to calculate the importance of the TF (default n=5)
#' @param buffer Distance buffer (in base pairs) for SNP-to-peak mapping (default = 500pb), which means each peak will be extended by 500bp upstream and 500bp downstream.
#' @return A list containing MORE_score (specificity*JSI) and MORE_perm_Pval (Monte Carlo P-values)
#' @export

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
  # Map peaks to genes
  peak2gene_strength <- peak2gene(grn_outputs)
  # Map SNPs to peaks using the specified distance buffer
  snp2peak_map <- snp2peak(snp_info, peak2gene_strength, buffer = buffer)

  #snp2peak_map <- snp2peak(snp_info_lymp, peak2gene_strength, buffer = buffer)
  #snp2peak_map <- snp2peak(snp_info_mono, peak2gene_strength, buffer = buffer)
  
  snp2peak_map$geneScores <- geneRiskScores$logP[match(snp2peak_map$Target,geneRiskScores$SYMBOL)]
  #sum(is.na(snp2peak_map$geneScores))
  
  # Calculate peak importance scores
  snp2peak_map <- getPeakScore(snp2peak_map)
  
  # Extract regulons containing SNP, peak, TF, target gene, and importance score
  regulons <- snp2peak_map[, c("snp_id", "peak_ids", "TF", "Target", "Importance_weighted")]
  
  # Step 2: TF list and cell types
  tf_list <- grn_outputs$tf_names  # List of transcription factors (TFs)
  all_celltype_names <- unique(target_scores[, "celltypes"])  # Unique cell types for analysis

  # Step 3: COSR and MORE analysis for each cell type and TF
  #Open progress bar
  pb <- txtProgressBar(style=3)
  start_time <- Sys.time() ##record start time

  total_run <- length(all_celltype_names)*length(tf_list)
  count <- 0
  
  # Initialize an empty data frame to store results
  # Initialize an empty data frame to store results
  all_regulon_results_df <- data.frame(
    RegulonID = character(),                # Unique identifier for each regulon
    RegulonName = character(),              # Name of the regulon
    SpecificityScore = numeric(),           # Specificity score
    SpecificityScore_p = numeric(),         # P-value for specificity score
    GeneRiskScore = numeric(),              # Gene risk score
    ImportanceWeightScore_p = numeric(),    # P-value for importance score
    RegulonScore = numeric(),               # Final regulon score
    RegulonScore_p = numeric(),             # P-value for regulon score
    Celltype = character(),                 # Cell type name
    stringsAsFactors = FALSE                # Prevent factors for character columns
  )
  
  #COSR and JSI interaction analysis
  for (i in seq_along(all_celltype_names)){

    for (j in seq_along(tf_list)){

      #extracting the gene list of each regulon
      
      #j=12
      #i=1
      Module <- regulons[which(regulons$TF == tf_list[j]),]
      Module_regulon <- c(unique(Module$TF),unique(Module$Target))

      #obtain the TF and target genes related importance_weighted scores in each regulon
      eachModule_Importance_score <- getRiskScore(Module,top_n = 5)

      #all specificty score and z score of all regulon genes
      target_scores_sub <- target_scores[which(target_scores[,3] == all_celltype_names[i]),]

      #extracting the specificity score of each regulon
      each_module_score <- target_scores_sub[!is.na(match(target_scores_sub[,1], Module_regulon)),]
      #each_module_score$anno <- rep("Gene",n_num,length(each_module_score[,1]))
      each_module_score$anno <- rep("Gene", length(each_module_score[,1]))
      each_module_score$anno[which(each_module_score[,1] == Module_regulon[1])] <- "TF"

      # add importance weighted score in each regulon
      each_module_score$Importance_weighted <- eachModule_Importance_score$Importance_weighted[match(each_module_score$genes,eachModule_Importance_score$Target)]


      # Trait-associated regulon score (TARS)
      TARS <- getRegulonScore(each_module_score,theta=0.5,alpha=0.5)
      
      # Run Monte Carlo (MC) Permutation analysis
      
      #tf_list <- tf_list[which(tf_list!=tf_list[j])] #removing the targeted TF as controls
      len_of_regulon <- length(Module_regulon)

      #get random importance weight, specificity, and risk scores for background genes
      #target_scores_background <- getRandomWeight(regulons,target_scores_sub)
      
      target_scores_background <- target_scores_sub
      #real_specificity <-  sample(target_scores_sub$scores[which(target_scores_sub$genes %in% regulons$Target)],length(regulons$snp_id),replace = T)   # Replace with your real specificity scores
      #real_specificity <-  sample(target_scores_sub$scores,length(regulons$snp_id),replace = T)   # Replace with your real specificity scores
      real_specificity <-  target_scores_sub$scores[which(target_scores_sub$genes %in% regulons$Target)]  # Replace with your real specificity scores
      
      real_importance <- regulons$Importance_weighted    # Replace with your real importance scores
      
      
      ##run MC permutation analysis
      perm_results <- replicate(1000,generate_random_regulon_scores(
                                                   tf_list,
                                                   target_scores_background,
                                                   real_specificity,
                                                   real_importance,
                                                   len_of_regulon,
                                                   theta=0.5,
                                                   alpha=0.5,
                                                   top_n = 5))
      
      # Extract scores from the matrix
      perm_specificity <- perm_results["SpecificityScore", ]
      perm_importance <- perm_results["ImportanceWeightScore", ]
      perm_regulon <- perm_results["RegulonScore", ]

      # transform numeric
      perm_specificity <- as.numeric(perm_specificity)
      perm_importance <- as.numeric(perm_importance)
      perm_regulon <- as.numeric(perm_regulon)
      
      # check and remove NA
      if (any(is.na(perm_specificity))) {
        warning("NA values found in perm_specificity. Removing them.")
        perm_specificity <- perm_specificity[!is.na(perm_specificity)]
      }
      
      if (any(is.na(perm_importance))) {
        warning("NA values found in perm_importance. Removing them.")
        perm_importance <- perm_importance[!is.na(perm_importance)]
      }
      
      if (any(is.na(perm_regulon))) {
        warning("NA values found in perm_regulon. Removing them.")
        perm_regulon <- perm_regulon[!is.na(perm_regulon)]
      }
      
      
      # Calculate p-values
      p_specificity <- (1 + sum(perm_specificity >= TARS$SpecificityScore)) / (1 + length(perm_specificity))
      p_importance <- (1 + sum(perm_importance >= TARS$GeneRiskScore)) / (1 + length(perm_importance))
      p_regulon <- (1 + sum(perm_regulon >= TARS$RegulonScore)) / (1 + length(perm_regulon))
      
      # Calculate z-score
      z_specificity <- (TARS$SpecificityScore - mean(perm_specificity)) / sd(perm_specificity)
      z_importance <- (TARS$GeneRiskScore - mean(perm_importance)) / sd(perm_importance)
      z_regulon <- (TARS$RegulonScore - mean(perm_regulon)) / sd(perm_regulon)
      
      
      # Collect the regulon results
      # Append the results to the data frame
      # Collect the regulon results
      all_regulon_results_df <- rbind(
        all_regulon_results_df, 
        data.frame(
          RegulonID = paste0("Regulon_", j),          # Unique identifier
          RegulonName = Module_regulon[1],           # Replace with the actual regulon name
          SpecificityScore = z_specificity,  # Specificity score
          SpecificityScore_p = p_specificity,        # P-value for specificity score
          GeneRiskScore = z_importance,        # Gene risk score
          ImportanceWeightScore_p = p_importance,    # P-value for importance score
          RegulonScore = z_regulon,          # Final regulon score
          RegulonScore_p = p_regulon,                # P-value for regulon score
          Celltype = all_celltype_names[i]           # Cell type name
        )
      )
      
      
      #Running
      print(paste("Running the regulon of ",tf_list[j], " for the cell type of ",all_celltype_names[i],sep = ""))

      #Running percent:
      count=count+1
      completed_percent <- count/total_run
      print(sprintf('Completed percent: %1.2f%%',100*completed_percent))
      #Real-time progress bar
      #print(paste("Runing percent: ",percent((i+j)/(length(all_celltype_names)*length(tf_list))),sep = ""))
      setTxtProgressBar(pb,(count)/total_run)

    }

  }
  

  # Format the numeric columns to display as decimals with a fixed number of digits (e.g., 4 digits)
  all_regulon_results_df$SpecificityScore <- format(all_regulon_results_df$SpecificityScore, nsmall = 4, scientific = FALSE)
  all_regulon_results_df$GeneRiskScore <- format(all_regulon_results_df$GeneRiskScore, nsmall = 4, scientific = FALSE)
  all_regulon_results_df$RegulonScore <- format(all_regulon_results_df$RegulonScore, nsmall = 4, scientific = FALSE)
  
  
  # define the significant regulons in each cell type
  all_regulon_results_df$Significance <- ifelse(
    all_regulon_results_df$SpecificityScore_p < 0.05 & 
      all_regulon_results_df$ImportanceWeightScore_p < 0.05 & 
      all_regulon_results_df$RegulonScore_p < 0.05,
    "Significant",
    "Nonsignificant"
  )
  
  
  #write.csv(all_regulon_results_df,file="all_regulon_results_df7_lymp_specificity7.csv",quote = F)
  #write.csv(all_regulon_results_df,file="all_regulon_results_df7_mono_specificity4.csv",quote = F)
  
  
  ##Record end time
  end_time <- Sys.time()

  #Close progress bar
  close(pb)

  #Calculating the running time
  print(paste("Running time:", round(as.numeric(difftime(end_time, start_time, units = "secs")), 2), "seconds"))

  #Outputs
  #out_results <- list(ctDRTF_score = Final_regulon_ct_score,
  #                     MC_p = Final_regulon_ct_mc_p,
  #                      regulon_specificity_s=Final_regulon_s)
  #
  MORE_results <- all_regulon_results_df

  #
  return(MORE_results)

}
