#' @title title: alternativeSpecificityScore
#' @description
#' This function calculates specificity scores for regulatory modules (regulons)
#' in single-cell RNA-seq data using one of the five alternative methods: AUCell, RSS, VAM, Seurat, or UCell.
#' @param single_cell A Seurat object containing single-cell RNA-seq data.
#' @param tf_list A vector of transcription factors (TFs) to evaluate.
#' @param regulons A data frame containing the regulatory relationships. It should have at least two columns:
#'             'TF' (transcription factor) and 'Target' (target genes).
#' @param alternative A string specifying the method to use for scoring. Options are:
#'                "AUCell", "RSS", "VAM", "Seurat", "UCell". Default is "AUCell".
#' @export
#'
alternativeSpecificityScore <- function(single_cell = single_cell,
                                        tf_list = tf_list,
                                        regulons = regulons,
                                        alternative = "AUCell") {

  # Initialize an empty data frame to store results
  regulon_score <- data.frame(celltypes = character(),   # Cell type identities
                              scores = numeric(),        # Regulatory scores
                              regulons = character(),    # Regulon (TF) names
                              stringsAsFactors = FALSE)

  # Loop through each transcription factor in tf_list
  for (j in seq_along(tf_list)) {
    # Extract regulon genes for the current TF
    Module <- regulons[regulons$TF == tf_list[j], ]      # Filter regulons by TF
    Module_regulon <- unique(c(Module$TF, Module$Target)) # Combine TF and its targets

    # Use the specified method to calculate scores
    if (alternative == "AUCell") {
      # AUCell method
      single_cell_matrix <- as.matrix(single_cell@assays$RNA@data)  # Convert RNA assay data to matrix
      cells_rankings <- AUCell_buildRankings(single_cell_matrix)   # Build rankings for cells
      cells_AUC <- AUCell_calcAUC(Module_regulon, cells_rankings,  # Calculate AUC for the regulon
                                  aucMaxRank = nrow(cells_rankings) * 0.05)
      AUCell_auc <- as.numeric(getAUC(cells_AUC)[1, ])             # Extract AUC scores
      auc_cell <- data.frame(celltypes = Idents(single_cell), scores = AUCell_auc) # Combine cell types with scores
      results_auc <- aggregate(scores ~ celltypes, auc_cell, mean) # Aggregate scores by cell type
      results_auc$regulons <- tf_list[j]                           # Add the current TF name
      regulon_score <- bind_rows(regulon_score, results_auc)       # Append results to the output

    } else if (alternative == "RSS") {
      # RSS method
      single_cell_matrix <- as.matrix(single_cell@assays$RNA@data)
      cells_rankings <- AUCell_buildRankings(single_cell_matrix)
      cells_AUC <- AUCell_calcAUC(Module_regulon, cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.05)
      AUCell_rss <- as.numeric(getRSS(cells_AUC, cell_anno = Idents(single_cell))[1, ])
      names(AUCell_rss) <- levels(Idents(single_cell))
      rss_cell <- data.frame(celltypes = names(AUCell_rss), scores = AUCell_rss, regulons = tf_list[j])
      regulon_score <- bind_rows(regulon_score, rss_cell)

    } else if (alternative == "VAM") {
      # VAM method
      library(SingleCellExperiment)
      library(qusage)
      sce_multi <- as.SingleCellExperiment(single_cell)            # Convert Seurat object to SCE
      geneset <- list(Module_regulon)                              # Create a gene set for the regulon
      names(geneset) <- tf_list[j]
      Module_regulon_gmt <- get_gmt(geneset)                       # Generate GMT format for the regulon
      gmt_file <- tempfile(fileext = ".gmt")                       # Create a temporary GMT file
      writeLines(Module_regulon_gmt, con = gmt_file)
      gs <- qusage::read.gmt(gmt_file)                             # Read the GMT file
      sce_multi <- importGeneSetsFromList(inSCE = sce_multi, geneSetList = gs, by = "rownames")
      sce_multi <- runVAM(inSCE = sce_multi, geneSetCollectionName = "GeneSetCollection", useAssay = "logcounts")
      VAM_results <- sce_multi@int_colData$reducedDims$VAM_GeneSetCollection_CDF[, 1] # Extract VAM scores
      vam_cell <- data.frame(celltypes = Idents(single_cell), scores = VAM_results)
      results_vam <- aggregate(scores ~ celltypes, vam_cell, mean)
      results_vam$regulons <- tf_list[j]
      regulon_score <- bind_rows(regulon_score, results_vam)

    } else if (alternative == "Seurat") {
      # Seurat AddModuleScore method
      DefaultAssay(single_cell) <- "RNA"
      single_cell <- AddModuleScore(single_cell, features = list(Module_regulon), name = "Module_regulon")
      add_cell <- data.frame(celltypes = Idents(single_cell), scores = single_cell$Module_regulon1)
      results_add <- aggregate(scores ~ celltypes, add_cell, mean)
      results_add$regulons <- tf_list[j]
      regulon_score <- bind_rows(regulon_score, results_add)

    } else if (alternative == "UCell") {
      # UCell method
      seu_matrix <- single_cell@assays$RNA@data
      UCell_scores <- ScoreSignatures_UCell(seu_matrix, features = list(Module_regulon))
      colnames(UCell_scores) <- "scores"
      ucell_cell <- data.frame(celltypes = Idents(single_cell), scores = UCell_scores)
      results_ucell <- aggregate(scores ~ celltypes, ucell_cell, mean)
      results_ucell$regulons <- tf_list[j]
      regulon_score <- bind_rows(regulon_score, results_ucell)
    }
  }

  # Remove rows with missing values
  regulon_score <- regulon_score[complete.cases(regulon_score), ]

  # Return the final result
  return(regulon_score)
}
