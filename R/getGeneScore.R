#'
#' @title Pre-process the gene-based genetic data
#'
#' @description
#' Extract and format gene-based genetic association results
#'
#' @param gene_info Data frame of gene-based genetic association results (MAGMA or FUMA format)
#' @importFrom dplyr %>% mutate arrange filter select
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @return A data frame with columns:
#'         - SYMBOL: Gene symbols
#'         - logP: -log10(P) values
#'         - ZSTAT: Z-scores
#' @export
#'
getGeneScore <- function(gene_info) {

  if (!is.data.frame(gene_info)) {
    stop("Input gene_info must be a data.frame or tibble.")
  }

  # Check if gene symbols are already provided
  if ("SYMBOL" %in% colnames(gene_info)) {
    # Process MAGMA-format data
    geneRiskScores <- gene_info %>%
      dplyr::mutate(logP = -log10(P)) %>%
      dplyr::arrange(desc(logP)) %>%
      dplyr::select(SYMBOL, logP, ZSTAT)
  } else {
    # Map gene identifiers to symbols using org.Hs.eg.db
    gene_info$SYMBOL <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = as.character(gene_info$GENE),
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )

    # Filter out rows with missing gene symbols and process the data
    geneRiskScores <- gene_info %>%
      dplyr::filter(!is.na(SYMBOL)) %>%
      dplyr::mutate(logP = -log10(P)) %>%
      dplyr::arrange(desc(logP)) %>%
      dplyr::select(SYMBOL, logP, ZSTAT)
  }

  return(geneRiskScores)
}
