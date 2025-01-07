#' @title Link genes to genomic peaks
#' @description Links genes to genomic peaks and computes Strength based on the GRN inference method.
#'
#' @param grn_outputs A list containing GRN output data, with a 'Regions' column specifying genomic peak ranges.
#' @param infer_method GRN inference method:
#'   - 'glm' (default): Generalized linear model
#'   - 'cv.glmnet': Regularized GLMs
#'   - 'brms': Bayesian regression models
#'   - 'bagging_ridge': Bagging ridge and Bayesian ridge (CellOracle)
#'   - 'xgb': XGBoost gradient-boosted Random Forest (SCENIC)
#'
#' @return A data frame mapping genes to genomic peaks with associated log-transformed p-values (Strength).
#'
#' @import dplyr
#' @import tidyr
#' @export

peak2gene <- function(grn_outputs, infer_method = 'glm') {
  # Extract GRN data
  grn_data <- grn_outputs$grn

  # Ensure the 'Regions' column is of character type
  if (!"Regions" %in% colnames(grn_data)) {
    stop("'Regions' column is missing in the GRN data")
  }
  grn_data$Regions <- as.character(grn_data$Regions)

  # Determine the column to use for Strength calculation
  strength_column <- switch(
    infer_method,
    glm = "Pval",
    cv.glmnet = "Corr",
    glmnet = "Corr",
    xgb = "Gain",
    stop("Invalid infer_method specified")
  )

  if (!strength_column %in% colnames(grn_data)) {
    stop(paste0("'", strength_column, "' column is missing in the GRN data"))
  }

  # Split the 'Regions' column by ";" and expand into multiple rows
  expanded_grn_data <- grn_data %>%
    tidyr::separate_rows(Regions, sep = ";") %>%
    dplyr::mutate(
      Strength = switch(
        infer_method,
        glm = -log10(Pval),
        cv.glmnet = (Corr + 1) / 2,
        glmnet = (Corr + 1) / 2,
        xgb = Gain + 0.00005
      )
    )

  return(expanded_grn_data)
}
