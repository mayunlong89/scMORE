#' @title Construct TF-gene regulatory network (GRN)
#' @description
#' Use single-cell multiomics data (scRNA-seq and scATAC-seq) to construct a global TF-peak-gene regulatory network.
#'
#' @param single_cell Input single-cell multiomics data (Seurat object)
#' @param n_targets Minimum number of targets required in a regulon (default = 5)
#' @param peak2gene_method Method for peak-to-gene linkage: 'Signac' or 'GREAT'
#' @param infer_method GRN inference method:
#'   - 'glm' (default): Generalized linear model
#'   - 'cv.glmnet': Regularized GLMs
#'   - 'brms': Bayesian regression models
#'   - 'bagging_ridge': Bagging ridge and Bayesian ridge (CellOracle)
#'   - 'xgb': XGBoost gradient-boosted Random Forest (SCENIC)
#'
#' @return A list containing:
#'   - grn: Filtered GRN with regulons meeting the gene threshold
#'   - tf_names: List of TFs with sufficient regulons
#'
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @import Pando
#' @export
createRegulon <- function(single_cell,
                          n_targets = 5,
                          peak2gene_method = 'Signac',
                          infer_method = 'glm',
                          tss_upstream = 100000,
                          tss_downstream=0,
                          exclude_exon_regions = TRUE,
                          conserved_regions = phastConsElements20Mammals.UCSC.hg38) {

  # Step 1: Select variable features from single-cell data
  single_cell <- Seurat::FindVariableFeatures(single_cell, assay = "RNA")
  data("phastConsElements20Mammals.UCSC.hg38", package = "scMORE")

  # Step 2: Initiate GRN object and select candidate regions
  single_cell <- Pando::initiate_grn(
    single_cell,
    peak_assay = "peaks",
    rna_assay = "RNA",
    exclude_exons = exclude_exon_regions,
    regions = conserved_regions
  )

  # Step 3: Scan candidate regions for TF binding motifs
  data("motifs", package = "scMORE")
  data("motif2tf", package = "scMORE")
  single_cell <- Pando::find_motifs(
    single_cell,
    pfm = motifs,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )

  # Step 4: Infer the gene regulatory network
  single_cell <- Pando::infer_grn(
    single_cell,
    peak_to_gene_method = peak2gene_method,
    method = infer_method,
    upstream=tss_upstream,
    downstream=tss_downstream,
    alpha = 0.5,
    family = 'gaussian',
    adjust_method = 'fdr',
    scale = FALSE,
    verbose = TRUE
  )

  # Step 5: Identify regulatory modules
  single_cell <- Pando::find_modules(
    single_cell,
    p_thresh = 0.1,
    nvar_thresh = 2,
    min_genes_per_module = 1,
    rsq_thresh = 0.05
  )

  # Step 6: Extract GRN modules (regulons)
  regulons <- Pando::NetworkModules(single_cell)

  # Step 7: Extract regulatory network data
  grn_data <- extract_grn(regulons, single_cell, infer_method)

  # Step 8: Filter regulons with at least `n_targets` target genes
  tf_gene_counts <- table(grn_data$TF)
  valid_tfs <- names(tf_gene_counts[tf_gene_counts >= n_targets])
  filtered_regulons <- grn_data[grn_data$TF %in% valid_tfs, ]

  # Step 9: Prepare outputs
  grn_outputs <- list(
    grn = filtered_regulons,
    tf_names = valid_tfs
  )

  return(grn_outputs)
}
