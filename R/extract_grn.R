#' @title Extract GRN Data
#' @description Helper function to extract GRN data based on inference method.
#'
extract_grn <- function(regulons, single_cell, infer_method) {
  if (infer_method %in% c("cv.glmnet", "glmnet")) {
    grn_data <- coef(single_cell)
    corr_info <- grn_data %>% dplyr::select(tf, target, region, corr)
    regulons_df <- as.data.frame(regulons@meta)
    regulons_with_corr <- regulons_df %>% dplyr::left_join(corr_info, by = c("tf", "target"))
    return(data.frame(
      TF = regulons_with_corr$tf,
      Target = regulons_with_corr$target,
      Regions = regulons_with_corr$regions,
      Corr = regulons_with_corr$corr
    ))
  } else if (infer_method == "xgb") {
    grn_data <- coef(single_cell)
    corr_info <- grn_data %>% dplyr::select(tf, target, region, gain)
    regulons_df <- as.data.frame(regulons@meta)
    regulons_with_corr <- regulons_df %>% dplyr::left_join(corr_info, by = c("tf", "target"))
    return(data.frame(
      TF = regulons_with_corr$tf,
      Target = regulons_with_corr$target,
      Regions = regulons_with_corr$regions,
      Gain = regulons_with_corr$gain
    ))
  } else if (infer_method == "glm") {
    return(data.frame(
      TF = regulons@meta$tf,
      Target = regulons@meta$target,
      Regions = regulons@meta$regions,
      Pval = regulons@meta$pval
    ))
  } else {
    stop("Invalid inference method")
  }
}
