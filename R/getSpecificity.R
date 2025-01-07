#' @title Calculate Cell-Type Specificity Scores of Nodes
#' @description
#' Calculate specificity scores (CTS) of TFs and target genes across all cell types.
#' @param single_cell Input single-cell data (Seurat object)
#' @param method Method used for calculating the cell type specificity of nodes (TFs or genes):
#'        - 'cosine' (default)
#'        - 'mean_specificity'
#'        - 'average'
#' @return A data frame with columns:
#'         - genes: Gene or TF names
#'         - scores: Specificity scores
#'         - celltypes: Cell type names
#' @export
getSpecificity <- function(single_cell, method = 'cosine') {
  # Step 1: Initialize variables
  num_genes <- nrow(single_cell) # Total number of genes in the dataset
  target_scores <- data.frame(genes = character(), scores = numeric(), celltypes = character(), stringsAsFactors = FALSE)

  # Step 2: Calculate specificity scores based on the selected method
  if (method == 'cosine') {
    # Using COSG for cosine similarity-based specificity
    celltype_markers <- COSG::cosg(
      single_cell,
      groups = "all",
      assay = "RNA",
      slot = "data",
      mu = 1,
      n_genes_user = num_genes
    )
    all_celltype_names <- colnames(celltype_markers$names) # Extract cell type names
    target_scores <- do.call(rbind, lapply(all_celltype_names, function(celltype) {
      data.frame(
        genes = celltype_markers$names[[celltype]],
        scores = celltype_markers$scores[[celltype]],
        celltypes = celltype,
        stringsAsFactors = FALSE
      )
    }))
  } else if (method == 'average') {
    all_celltype_names <- as.vector(unique(Idents(single_cell))) # Extract cell type names
    celltype_averages <- as.data.frame(AverageExpression(single_cell)$RNA) # Average expression values
    for (celltype in all_celltype_names) {
      temp <- data.frame(
        genes = rownames(celltype_averages),
        scores = celltype_averages[[celltype]],
        celltypes = celltype,
        stringsAsFactors = FALSE
      )
      target_scores <- rbind(target_scores, temp)
    }
  } else if (method == "mean_specificity") {
    # Step 1: Calculate the average expression of each gene for each cell type
    celltype_averages <- as.data.frame(AverageExpression(single_cell)$RNA)

    # Step 2: Add the average expression of each gene across all cell types
    celltype_averages$Overall_Mean <- rowMeans(celltype_averages, na.rm = TRUE) # Handle NA values
    celltype_averages$Overall_Mean[celltype_averages$Overall_Mean == 0 | is.na(celltype_averages$Overall_Mean)] <- 1e-6

    # Step 3: Calculate specificity for each cell type
    specificity <- sweep(celltype_averages[, -ncol(celltype_averages)], 1, celltype_averages$Overall_Mean, "/")
    # Step 4: Convert specificity matrix to a data frame
    specificity_df <- as.data.frame(specificity)
    # Step 5: Add gene names as a new column
    specificity_df$genes <- rownames(specificity_df)

    # Step 6: Reshape the data frame to a long format
    target_scores <- do.call(rbind, lapply(names(specificity_df)[-ncol(specificity_df)], function(celltype) {
      data.frame(
        genes = specificity_df$genes,
        celltypes = celltype,
        scores = specificity_df[[celltype]],
        stringsAsFactors = FALSE
      )
    }))
  } else {
    stop("Invalid specificity_method. Choose 'cosine', 'average', or 'mean_specificity'.")
  }

  # Step 3: Replace specificity scores of -1 with 0
  target_scores$scores[target_scores$scores == -1] <- 0

  # Step 4: Log10-transformation and max-min scale the specificity scores
  target_scores$scores <- log10(target_scores$scores + 1e-6)
  target_scores$scores <- (target_scores$scores - min(target_scores$scores, na.rm = TRUE)) /
    (max(target_scores$scores, na.rm = TRUE) - min(target_scores$scores, na.rm = TRUE))

  # Step 5: Return the final data frame
  return(target_scores)
}
