#' @title Calculate Cell-Type Specificity Scores of Nodes
#'
#' @description
#' Calculate specificity scores (CTS) of TFs and target genes across all cell types.
#' @param single_cell Input single-cell data (Seurat object)
#' @param specificity_method Method used for calculating the cell type specificity of nodes (TFs or genes):
#'        'cosine' (default) or 'average'.
#' @return A data frame with columns:
#'         - genes: Gene or TF names
#'         - scores: Specificity scores
#'         - celltypes: Cell type names
#' @export
#'
getSpecificity <- function(single_cell, specificity_method = 'cosine') {
  # Step 1: Initialize variables
  num_genes <- nrow(single_cell) # Total number of genes in the dataset
  target_scores <- data.frame(genes = character(), scores = numeric(), celltypes = character(), stringsAsFactors = FALSE)

  # Step 2: Calculate specificity scores based on the selected method
  if (specificity_method == 'cosine') {
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

    # Combine specificity scores for all cell types into a single data frame
    target_scores <- do.call(rbind, lapply(all_celltype_names, function(celltype) {
      data.frame(
        genes = celltype_markers$names[[celltype]],
        scores = celltype_markers$scores[[celltype]],
        celltypes = celltype,
        stringsAsFactors = FALSE
      )
    }))

  } else if (specificity_method == 'average') {
    # Calculate average-based specificity scores
    all_celltype_names <- as.vector(unique(Idents(single_cell))) # Extract cell type names
    celltype_averages <- AverageExpression(single_cell)$RNA # Average expression values

    for (celltype in all_celltype_names) {
      temp <- data.frame(
        genes = rownames(celltype_averages),
        scores = celltype_averages[[celltype]],
        celltypes = celltype,
        stringsAsFactors = FALSE
      )
      target_scores <- rbind(target_scores, temp)
    }
  } else {
    stop("Invalid specificity_method. Choose 'cosine' or 'average'.")
  }

  # Step 3: Replace specificity scores of -1 with 0
  target_scores$scores[target_scores$scores == -1] <- 0

  # Step 4: Log10-transformation and max-min scale the specificity scores
  target_scores$scores <- log10(target_scores$scores + 1e-6)
  target_scores$scores <- (target_scores$scores - min(target_scores$scores)) / (max(target_scores$scores) - min(target_scores$scores))

  # Step 5: Return the final data frame
  return(target_scores)
}
