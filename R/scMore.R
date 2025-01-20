#' @title Single-Cell MultiOmics Regulon Enrichment (scMORE)
#'
#' @description
#' Main function for identifying disease-relevant regulons using single-cell data.
#'
#' @param single_cell The input single-cell data, typically a Seurat or SingleCellExperiment object.
#' @param snp_info GWAS summary statistics. Must include the following columns:
#'                 `chr` (chromosome), `pos` (position), `rsID` (SNP ID), `beta` (effect size),
#'                 and `p-value` (significance level).
#' @param n_targets Minimum number of target genes required in a given regulon. Default: 5.
#' @param gene_info Gene-based genetic association results from MAGMA or FUMA.
#'                  The function supports results in standard MAGMA/FUMA formats.
#' @param method Method used for calculating the cell type specificity of nodes (TFs or genes):
#'        - 'cosine' (default)
#'        - 'mean_specificity'
#'        - 'average'
#' @param perm_n Number of Monte Carlo permutations to perform for significance testing. Default: 1000.
#' @param theta Weighting factor to adjust the specificity score of target genes within each regulon.
#'              Range: 0.1 to 1. Default: 0.5.
#' @param alpha Flexibility parameter for penalization in the inference model. Default: 1.
#' @param buffer Numeric value specifying the flanking region size for genomic peaks.
#'               Extends the peak range upstream and downstream. For example, `buffer = 500`
#'               extends 500 bp on both sides of a peak. Default: 500 bp.
#' @param p1 Threshold for statistical significance of the cell type-specificity score
#'           for each regulon. Default: 0.05.
#' @param p2 Threshold for statistical significance of the genetic risk score
#'           for each regulon. Default: 0.05.
#' @param p3 Threshold for statistical significance of the trait-associated regulon score (TRS).
#'           Default: 0.05.
#' @param peak2gene_method Method for mapping peaks to genes. Options: `'Signac'` or `'GREAT'`.
#' @param infer_method Method for inference modeling. Options:
#'                     - Generalized Linear Model (`'glm'`)
#'                     - Regularized GLM (`'glmnet'`, `'cv.glmnet'`)
#'                     - Bayesian regression models (`'brms'`)
#'                     - Bagging ridge and Bayesian ridge models (`'bagging_ridge'`,CellOracle)
#'                     - XGBoost gradient-boosted Random Forest (`'xgb'`, SCENIC)
#' @param top_n Number of top targets for each TF used to calculate its importance. Default: 5.
#' @param nSeed Random seed for reproducibility (i.e., `set.seed()`). Default: 1234.
#'
#' @return A list containing:
#'         - Disease-relevant regulons
#'         - Their genetic and cell-type associations
#'
#' @export
#'
scMore <- function(single_cell,
                   snp_info,
                   gene_info,
                   n_targets = 5,
                   perm_n = 1000,
                   theta = 0.5,
                   alpha = 1,
                   buffer = 500,
                   top_n = 5,
                   p1 = 0.05,
                   p2 = 0.05,
                   p3 = 0.05,
                   peak2gene_method = 'Signac',
                   infer_method = 'glm',
                   method = 'cosine',
                   nSeed = 1234) {

  # Set random seed for reproducibility
  set.seed(nSeed)

  # Step 1: Construct the global TF-gene regulatory network
  # Identify transcription factors (TFs) and their target genes
  message("Step 1: Constructing the global TF-gene regulatory network...")
  grn_outputs <- createRegulon(single_cell,
                               n_targets = n_targets,
                               peak2gene_method = peak2gene_method,
                               infer_method = infer_method)

  # Step 2: Calculate specificity scores for target genes
  # Evaluate the relevance of target genes to specific regulons
  message("Step 2: Calculating specificity scores for target genes...")
  target_scores <- suppressWarnings(getSpecificity(single_cell,
                                                   method = method))

  # Step 3: Get gene-level association scores
  message("Step 3: Retrieving gene-level association scores...")
  geneRiskScores <- getGeneScore(gene_info)

  # Step 4: Identify cell type-specific regulons relevant to disease
  # Integrate GRN outputs, target specificity scores, and gene association data
  message("Step 4: Identifying cell type-specific regulons relevant to disease...")
  regulon2disease_results <- regulon2disease(
    grn_outputs = grn_outputs,
    target_scores = target_scores,
    geneRiskScores = geneRiskScores,
    snp_info = snp_info,
    perm_n = perm_n,
    theta = theta,
    alpha = alpha,
    buffer = buffer,
    infer_method = infer_method,
    top_n = top_n
  )

  # Step 5: Collect all relevant results
  message("Step 5: Collecting all relevant results...")
  scMore_results <- list(
    scMore_grn = grn_outputs,
    scMore_celltype_specificity = geneRiskScores,
    scMore_trait_results = regulon2disease_results
  )

  # Step 6: Return final results
  message("Analysis completed successfully.")
  return(scMore_results)
}
