#'
#' @title Single-Cell MultiOmics REgulon discovery (scMORE)
#'
#' @description
#' Main function for identifying disease-relevant regulons using single-cell data
#'
#' @param single_cell The input single-cell data
#' @param snp_info GWAS summary statistics, must include columns: chr, pos, rsID, beta, p-value
#' @param n_targets Minimum number of targets required in a given regulon (default = 5)
#' @param gene_info Gene-based genetic association results based on MAGMA or FUMA
#'                  scMORE support the format of results from MAGMA and FUMA
#' @param top_genes Number of top-ranked disease-specific genes for JSI enrichment (default = 500)
#' @param perm_n Number of Monte Carlo permutations for significance testing (default = 1000)
#' @param theta Weight to adjust the specificity score of target genes in each regulon
#'              Range: 0.1 ~ 1 (default = 0.5)
#' @param pow Soft-power factor to emphasize strong specificity (default = 1)
#' @param mo Statistical model selection:
#'           - mo = 1: Genetic weight model
#'           - mo = 0: Non-genetic weight model
#'@param buffer Numeric value specifying the flanking sizes (default = 0)
#'          - extend peak range (upstream and downstream + 50bp)
#'          - buffer = 50bp, represent add 50 bp to both upstream and downstream of peak
#' @param nSeed set the seed for random sampling, i.e., set.seed()
#'
#' @return A list containing disease-relevant regulons and their associations
#' @export
#'
scMore <- function(single_cell,
                    snp_info,
                    gene_info,
                    n_targets= 5,
                    top_genes = 500,
                    perm_n = 1000,
                    theta=0.5,
                    pow=1,
                    mo=1,
                    buffer=50,
                    nSeed=1234){

  # Set random seed for reproducibility
  set.seed(nSeed)

  # Step 1: Construct the global TF-gene regulatory network
  # This step identifies transcription factors (TFs) and their target genes
  grn_outputs <- createRegulon(single_cell, n_targets)

  # Step 2: Calculate specificity scores for target genes
  # This step evaluates the relevance of target genes to specific regulons

  target_scores <- getSpecificity(single_cell)

  #get gene-level association scores
  geneRiskScores <- getGeneScore(gene_info)

  # Step 3: Identify cell type-specific regulons relevant to disease
  # This step integrates to determine disease associations

  regulon2disease_results <- regulon2disease(grn_outputs,
                                   target_scores,
                                   snp_info,
                                   geneRiskScores)

  scMore_results <- list(scMore_grn=grn_outputs,
                         scMore_celltype_specificity=geneRiskScores,
                         scMore_trait_results=regulon2disease_results)

  # Step 4: Return final results
  return(scMore_results)
}
