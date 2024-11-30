#'
#' @title Pre-process the gene-based genetic data
#'
#' @description
#' Extract and format gene-based genetic association results
#'
#' @param gene_info Data frame of gene-based genetic association results (MAGMA or FUMA format)
#' @return A data frame with columns:
#'         - SYMBOL: Gene symbols
#'         - logP: -log10(P) values
#'         - ZSTAT: Z-scores
#' @export
#'
getGeneScore <- function(gene_info) {
  # Check if gene symbols are already provided
  if ("SYMBOL" %in% colnames(gene_info)) {
    # Process MAGMA-format data
    geneRiskScores <- gene_info %>%
      mutate(logP = -log10(P)) %>%
      arrange(desc(logP)) %>%
      dplyr::select(SYMBOL, logP, ZSTAT)
  } else {
    # Map gene identifiers to symbols using org.Hs.eg.db
    gene_info$SYMBOL <- mapIds(
      org.Hs.eg.db,
      keys = as.character(gene_info$GENE),
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )

    # Filter out rows with missing gene symbols and process the data
    geneRiskScores <- gene_info %>%
      filter(!is.na(SYMBOL)) %>%
      mutate(logP = -log10(P)) %>%
      arrange(desc(logP)) %>%
      dplyr::select(SYMBOL, logP, ZSTAT)
  }

  return(geneRiskScores)
}
