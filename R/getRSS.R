#' @title getRSS
#' @description Calculates the Regulon Specificity Score (RSS).
#' The RSS quantifies the specificity of a regulon for a specific cell type based on AUCell results.
#'
#' @param cells_AUC A matrix of normalized AUCell results, where rows represent regulons and columns represent cells.
#' @param cell_anno A vector annotating each cell with its corresponding cell type.
#' @param cellTypes A vector of cell types to calculate RSS for. If NULL, it defaults to the unique values in `cell_anno`.
#'
#' @return A matrix where rows are regulons and columns are cell types, with the RSS values.
#'
#' @seealso
#' The RSS was first used by Suo et al. in:
#' "Revealing the Critical Regulators of Cell Identity in the Mouse Cell Atlas."
#' Cell Reports (2018). doi: 10.1016/j.celrep.2018.10.045.
#'
#' @examples
#' # Example usage:
#' # rss <- getRSS(cells_AUC, cell_anno)
#'
#' @export
getRSS <- function(cells_AUC, cell_anno, cellTypes = NULL) {
  # Check for NA values in cell annotations
  if (any(is.na(cell_anno))) stop("NAs in annotation")

  # If input is auCellResults, extract the AUC matrix
  if (any(class(cells_AUC) == "aucellResults")) cells_AUC <- getAUC(cells_AUC)

  # Normalize AUCell scores
  normAUC <- cells_AUC / rowSums(cells_AUC)

  # Set cell types to unique values in annotation if not provided
  if (is.null(cellTypes)) cellTypes <- unique(cell_anno)

  # Parallelized calculation of RSS if BiocParallel is available
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    ctapply <- BiocParallel::bplapply
  } else {
    ctapply <- lapply
  }

  # Calculate RSS for each cell type and regulon
  rss <- ctapply(cellTypes, function(thisType) {
    sapply(rownames(normAUC), function(thisRegulon) {
      pRegulon <- normAUC[thisRegulon, ]
      pCellType <- as.numeric(cell_anno == thisType)
      pCellType <- pCellType / sum(pCellType)
      .calcRSS.oneRegulon(pRegulon, pCellType)
    })
  })

  # Combine results into a matrix
  rss <- do.call(cbind, rss)
  colnames(rss) <- cellTypes
  return(rss)
}
