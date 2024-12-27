#' @title Single-Cell MultiOmics REgulon discovery (scMORE)
#'
#' @description
#' Main function for identifying disease-relevant regulons using single-cell data.
#'
#' @param single_cell The input single-cell data.
#' @param snp_info GWAS summary statistics, must include columns: chr, pos, rsID, beta, p-value.
#' @param n_targets Minimum number of targets required in a given regulon (default = 5).
#' @param gene_info Gene-based genetic association results based on MAGMA or FUMA,
#'                  (scMORE supports the format of results from MAGMA and FUMA).
#' @param perm_n Number of Monte Carlo permutations for significance testing (default = 1000).
#' @param theta Weight to adjust the specificity score of target genes in each regulon,
#'              Range: 0.1 ~ 1 (default = 0.5).
#' @param alpha Flexibility in penalization (default = 1).
#' @param buffer Numeric value specifying the flanking sizes (default = 500bp).
#'          - Extend peak range (upstream and downstream by 500bp)
#'          - buffer = 500bp represents adding 500bp to both upstream and downstream of a peak
#' @param peak2gene_method Two methods to choose: 'Signac' and 'GREAT'.
#' @param infer_method Users can choose different inference methods: GLMs ('glm'),
#'                    regularized GLMs('glmnet','cv.glmnet'), or Bayesian regression models('brms').
#' @param top_n Top n targets of each TF to calculate the importance of the TF (default = 5).
#' @param nSeed Set the seed for random sampling, i.e., set.seed() (default = 1234).
#'
#' @return A list containing disease-relevant regulons and their associations.
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
                   nSeed = 1234) {

  # Set random seed for reproducibility
  set.seed(nSeed)

  # Step 1: Construct the global TF-gene regulatory network
  # Identify transcription factors (TFs) and their target genes
  message("Step 1: Constructing the global TF-gene regulatory network...")
  grn_outputs <- createRegulon(single_cell, n_targets = n_targets)

  # Step 2: Calculate specificity scores for target genes
  # Evaluate the relevance of target genes to specific regulons
  message("Step 2: Calculating specificity scores for target genes...")
  target_scores <- suppressWarnings(getSpecificity(single_cell))

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
