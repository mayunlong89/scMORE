#'
#' @title get specificity in each cell type
#'
#' @description
#' Calculate specificity scores of TFs and target genes across all cell types
#'
#' @param single_cell Input single-cell data (Seurat object)
#' @return A data frame with columns:
#'         - genes: Gene or TF names
#'         - scores: Specificity scores
#'         - celltypes: Cell type names
#' @export
#'
getSpecificity <- function(single_cell) {

  # Step 1: Calculate specificity scores using COSG
  num_genes <- nrow(single_cell)  # Total number of genes in the dataset

  celltype_markers <- COSG::cosg(
    single_cell,
    groups = "all",
    assay = "RNA",
    slot = "data",
    mu = 1,
    n_genes_user = num_genes
  )

  # Step 2: Extract cell type names
  all_celltype_names <- colnames(celltype_markers$names)

  # Step 3: Combine specificity scores for all cell types into a single data frame
  target_scores <- do.call(rbind, lapply(all_celltype_names, function(celltype) {
    data.frame(
      genes = celltype_markers$names[[celltype]],
      scores = celltype_markers$scores[[celltype]],
      celltypes = celltype,
      stringsAsFactors = FALSE
    )
  }))

  # Step 4: Replace specificity scores of -1 with 0
  target_scores$scores[target_scores$scores == -1] <- 0
  
  # log10-transformation and max-min scale the target specificity score
  target_scores$scores <- log10(target_scores$scores + 1e-6)
  target_scores$scores <- max_min_scale(target_scores$scores) 
  
  # Step 5: Return the final data frame
  return(target_scores)
}
