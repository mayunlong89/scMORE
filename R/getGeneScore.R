#'
#' @title Pre-process the gene-based genetic data
#'
#' @description
#' Extract and format gene-based genetic association results
#'
#' @param gene_info Data frame of gene-based genetic association results (MAGMA or FUMA format)
#' @importFrom dplyr %>% mutate arrange filter select
#' @import dplyr
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @return A data frame with columns:
#'         - SYMBOL: Gene symbols
#'         - logP: -log10(P) values
#'         - ZSTAT: Z-scores
#' @export
#'
getGeneScore <- function(gene_info, p_floor = 1e-300) {
  if (!is.data.frame(gene_info))
    stop("Input gene_info must be a data.frame or tibble.")

  process_df <- function(df) {
    df %>%
      # 1. 先把 P=0 或极小值替换为 p_floor
      dplyr::mutate(P_adj = ifelse(P <= 0 | is.na(P), p_floor, P),
                    logP  = -log10(P_adj)) %>%
      dplyr::arrange(dplyr::desc(logP)) %>%
      dplyr::select(SYMBOL, logP, ZSTAT)
  }

  if ("SYMBOL" %in% colnames(gene_info)) {
    geneRiskScores <- process_df(gene_info)
  } else {
    gene_info$SYMBOL <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys      = as.character(gene_info$GENE),
      column    = "SYMBOL",
      keytype   = "ENTREZID",
      multiVals = "first"
    )
    geneRiskScores <- gene_info %>%
      dplyr::filter(!is.na(SYMBOL)) %>%
      process_df()
  }

  return(geneRiskScores)
}
