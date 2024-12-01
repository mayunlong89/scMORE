#' @title Identify cell type-specific regulons relevant to disease using COSR and MORE methods
#'
#' @param grn_outputs GRN outputs containing TF-gene relationships
#' @param target_scores Matrix of genes and TF specificity scores across cell types (columns: genes, scores, celltypes)
#' @param snp_info Information about SNPs for the analysis
#' @param geneRiskScores MAGMA-based gene association results (columns: SYMBOL, logP, ZSTAT)
#' @param perm_n Number of permutations for Monte Carlo simulation (default = 1000)
#' @param top_genes Number of top-ranked genes used for JSI calculation (default = 500)
#' @param theta Weight for integrating TF and gene scores (range: 0.1~1, default = 0.5)
#' @param pow Power parameter to emphasize specificity difference (default = 1)
#' @param mo Statistical model selection: mo = 1 for genetic weight model, mo = 0 for non-genetic model
#' @param buffer Distance buffer (in base pairs) for SNP-to-peak mapping (default = 50)
#' @return A list containing MORE_score (specificity*JSI) and MORE_perm_Pval (Monte Carlo P-values)
#' @export

regulon2disease <- function(grn_outputs,
                            target_scores,
                            snp_info,
                            geneRiskScores,
                            perm_n = 1000,
                            top_genes = 500,
                            theta = 0.5,
                            pow = 1,
                            mo = 1,
                            buffer = 50) {

  # Step 1: Map SNPs to TF-peaks and target genes
  # Map peaks to genes
  peak2gene_strength <- peak2gene(grn_outputs)
  # Map SNPs to peaks using the specified distance buffer
  snp2peak_map <- snp2peak(snp_info, peak2gene_strength, buffer = buffer)
  # Calculate peak importance scores
  snp2peak_map <- getPeakScore(snp2peak_map)

  # Extract regulons containing SNP, peak, TF, target gene, and importance score
  regulons <- snp2peak_map[, c("snp_id", "peak_ids", "TF", "Target", "Importance_weighted")]

  # Step 2: Initialize final result storage
  tf_list <- grn_outputs$tf_names  # List of transcription factors (TFs)
  Final_regulon_score <- data.frame(ID_regulon = tf_list)  # Store regulon scores
  Final_regulon_MORE_score <- data.frame(ID_regulon = tf_list)  # Store MORE scores (specificity * JSI)
  Final_regulon_MORE_perm_p <- data.frame(ID_regulon = tf_list)  # Store Monte Carlo permutation P-values
  all_celltype_names <- unique(target_scores[, "celltypes"])  # Unique cell types for analysis

  # Normalize target scores within each cell type
  for (m in all_celltype_names) {
    idx <- which(target_scores$celltypes == m)
    target_scores$scores[idx] <- max_min_scale(target_scores$scores[idx])  # Apply min-max scaling
  }

  # Step 3: COSR and MORE analysis for each cell type and TF
  #Open progress bar
  pb <- txtProgressBar(style=3)
  start_time <- Sys.time() ##record start time

  total_run <- length(all_celltype_names)*length(tf_list)
  count <- 0

  #COSR and JSI interaction analysis
  for (i in seq_along(all_celltype_names)){
    regulon_MORE_score <- c()
    regulon_MORE_perm_p <- c()
    regulon_score_all <-c()

    for (j in seq_along(tf_list)){

      #extracting the gene list of each regulon

      Module <- regulons[which(regulons$TF == tf_list[j]),]
      Module_regulon <- c(unique(Module$TF),unique(Module$Target))

      #obtain the TF and target genes related importance_weighted scores in each regulon
      eachModule_Importance_score <- getModuleScore(Module)

      #all specificty score and z score of all regulon genes
      target_scores_sub <- target_scores[which(target_scores[,3] == all_celltype_names[i]),]

      #annotating magma z-score
      geneRiskScores_sub <- geneRiskScores[which(geneRiskScores$SYMBOL %in% target_scores_sub$genes),]
      geneRiskScores_sub <- geneRiskScores_sub[!duplicated(geneRiskScores_sub[,c("SYMBOL")]),]

      #overlap target_scores regulon genes with magma genes
      target_scores_sub <- target_scores_sub[which(target_scores_sub$genes %in% geneRiskScores_sub$SYMBOL),]

      #match() function
      #data_set for all regulon genes specificity and z scores
      target_scores_sub$magma_zscore <- geneRiskScores_sub$ZSTAT[match(target_scores_sub$genes, geneRiskScores_sub$SYMBOL)]


      #extracting the specificity score of each regulon
      each_module_score <- target_scores_sub[!is.na(match(target_scores_sub[,1], Module_regulon)),]
      #each_module_score$anno <- rep("Gene",n_num,length(each_module_score[,1]))
      each_module_score$anno <- rep("Gene", length(each_module_score[,1]))
      each_module_score$anno[which(each_module_score[,1] == Module_regulon[1])] <- "TF"

      # add importance weighted score in each regulon
      each_module_score$Importance_weighted <- eachModule_Importance_score$Importance_weighted[match(each_module_score$genes,eachModule_Importance_score$Target)]


      #annotation MAGMA z-score
      #each_module_score<- each_module_score[which(each_module_score$genes %in%geneRiskScores$SYMBOL),]
      #each_module_score$magma_zscore <-geneRiskScores$ZSTAT[which(geneRiskScores$SYMBOL %in% each_module_score$genes)]


      #Calculating the module specificity score for TF in each regulon
      #tf_s_z <- each_module_score[,c("adj_score","magma_zscore")][which(each_module_score[,1] == Module_regulon[1]),]
      tf_Score_Zscore <- each_module_score[,c("scores","magma_zscore","Importance_weighted")][which(each_module_score[,1] == Module_regulon[1]),]
      tf_combined_score <- as.numeric((tf_Score_Zscore[1])^pow*(tf_Score_Zscore[2]*tf_Score_Zscore[3])^mo)

      if(is.na(tf_combined_score)){

        gene_Score_Zscore <- each_module_score[,c("scores","magma_zscore","Importance_weighted")][which(each_module_score[,1] != Module_regulon[1]),]
        gene_combined_score <- as.numeric((gene_Score_Zscore[,1])^pow*(gene_Score_Zscore[,2]*gene_Score_Zscore[,3])^mo)
        average_score <- sum(gene_combined_score)/(length(gene_combined_score)+1)

        tf_combined_score <- 0

        #theta = 0.5  #theta range from 0.1 ~ 1, default set to 0.5
        regulon_score <- as.numeric(tf_combined_score) + as.numeric(theta*average_score) #regulon-specific score for each cell type


      } else{

        #Calculating the module specificity score for genes in each regulon
        #gene_s_z <- each_module_score[,c("adj_score","magma_zscore")][which(each_module_score[,1] != Module_regulon[1]),]
        gene_Score_Zscore <- each_module_score[,c("scores","magma_zscore","Importance_weighted")][which(each_module_score[,1] != Module_regulon[1]),]
        gene_combined_score <- as.numeric((gene_Score_Zscore[,1])^pow*(gene_Score_Zscore[,2]* gene_Score_Zscore[,3])^mo)
        average_score <- mean(gene_combined_score)

        #theta = 0.5  #theta range from 0.1 ~ 1, default set to 0.5
        regulon_score <- as.numeric(tf_combined_score) + as.numeric(theta*average_score) #regulon-specific score for each cell type

      }

      #Sum
      regulon_score_all <- c(regulon_score_all,regulon_score)

      #Calculating the Jaccard Similarity Index (JSI)
      top_ranked_genes <-geneRiskScores$SYMBOL[1:top_genes]
      inter_genes_n <- length(intersect(top_ranked_genes,Module_regulon))
      union_genes_n <- length(union(top_ranked_genes,Module_regulon))
      JSI_score <- (inter_genes_n+1)/(union_genes_n+1) # Jaccard similarity index

      #Interaction: specificity*JSI for each regulon-disease link
      MORE_score <- regulon_score*JSI_score
      #MORE_score <- regulon_score

      #print(paste0("Regulon ",tf_list[j]," ctDRTF score is: ", MORE_score, sep=""))


      #Monte Carlo permutation for random specificity*JSI scores
      #Function: permutation()
      #perm_n = 1000
      tf_list_1 <- tf_list[which(tf_list!=tf_list[j])] #removing the targeted TF as controls
      len_of_regulon <- length(Module_regulon)
      all_genes <- geneRiskScores$SYMBOL[!is.na(geneRiskScores$SYMBOL)]

      #get random importance weight, specificity, and risk scores for background genes
      #target_scores_background <- getRandomWeight(regulons,target_scores_sub)

      random_scores <- getRandomScore(target_scores_sub)
      target_scores_background <- getRandomWeight(regulons,random_scores)



      ##run permutation analysis
      perm_results <- replicate(perm_n,permutation(target_scores_background,
                                                   tf_list_1,
                                                   len_of_regulon,
                                                   all_genes,
                                                   top_genes,
                                                   theta,
                                                   pow,
                                                   mo))


      dat <- as.numeric(perm_results)
      hist(dat, breaks = 50, col = "skyblue", main = "Improved Distribution")
      abline(v=MORE_score,col="red")

      #Calculating the MC p-values
      perm_p <- (1+length(perm_results[perm_results> MORE_score]))/(1+length(perm_results))


      #Running
      print(paste("Running the regulon of ",tf_list[j], " for the cell type of ",all_celltype_names[i],sep = ""))

      #Running percent:
      count=count+1
      completed_percent <- count/total_run
      print(sprintf('Completed percent: %1.2f%%',100*completed_percent))
      #Real-time progress bar
      #print(paste("Runing percent: ",percent((i+j)/(length(all_celltype_names)*length(tf_list))),sep = ""))
      setTxtProgressBar(pb,(count)/total_run)


      #Saving results
      regulon_MORE_score <- c(regulon_MORE_score,MORE_score)
      regulon_MORE_perm_p <- c(regulon_MORE_perm_p,perm_p)


    }

    #Collecting specificity scores
    regulon_score_all<- as.data.frame(regulon_score_all)
    names(regulon_score_all) <- all_celltype_names[i]
    Final_regulon_score <- cbind(Final_regulon_score,regulon_score_all)

    #Collecting MC P values
    regulon_MORE_perm_p<- as.data.frame(regulon_MORE_perm_p)
    names(regulon_MORE_perm_p) <- all_celltype_names[i]
    Final_regulon_MORE_perm_p <- cbind(Final_regulon_MORE_perm_p,regulon_MORE_perm_p)

    #Normalization

    regulon_MORE_score_norm <- (regulon_MORE_score - mean(regulon_MORE_score))/sd(regulon_MORE_score)
    #regulon_MORE_score_norm <- regulon_MORE_score

    #Alternative normalized method
    #max-min normalization
    #regulon_ct_score_norm <- max_min_scale(regulon_ct_score)

    #Collecting specificity*JSI for each regulon-disease link
    regulon_MORE_score_norm <- as.data.frame(regulon_MORE_score_norm)
    names(regulon_MORE_score_norm) <- all_celltype_names[i]
    Final_regulon_MORE_score <- cbind(Final_regulon_MORE_score,regulon_MORE_score_norm)

  }

  ##Record end time
  end_time <- Sys.time()

  #Close progress bar
  close(pb)

  #Calculating the running time
  end_time <- Sys.time()
  print(paste("Running time:", round(as.numeric(difftime(end_time, start_time, units = "secs")), 2), "seconds"))

  #Outputs
  #out_results <- list(ctDRTF_score = Final_regulon_ct_score,
  #                     MC_p = Final_regulon_ct_mc_p,
  #                      regulon_specificity_s=Final_regulon_s)
  #
  MORE_results <- list(MORE_score = Final_regulon_MORE_score,
                       MORE_perm_Pval = Final_regulon_MORE_perm_p)

  #
  return(MORE_results)

}
