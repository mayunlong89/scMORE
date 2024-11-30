#'
#' @title Construct TF-gene regulatory network (GRN)
#'
#' @description
#' Using single-cell multiomics data (scRNA-seq and scATAC-seq data) to construct a global TF-peak-gene regulatory network
#'
#'
#' @param single_cell Input single-cell multiomics data (Seurat object)
#' @param n_targets Minimum number of targets required in a regulon (default = 5)
#' @param peak2gene_method Two methods to choose: 'Signac' and 'GREAT'
#' @param infer_method You can choose different inference methods: GLMs ('glm'),
#'                    regularized GLMs('glmnet','cv.glmnet'), or Bayesian regression models('brms')
#' @return A list containing:
#'         - grn: Filtered GRN with regulons meeting the gene threshold
#'         - tf_names: List of TFs with sufficient regulons
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @export
#'
createRegulon <- function(single_cell,
                          n_targets = 5,
                          peak2gene_method='Signac',
                          infer_method='glm') {

  # Step 1: Select variable features from single-cell data
  single_cell <- Seurat::FindVariableFeatures(single_cell, assay = "RNA")


  data("phastConsElements20Mammals.UCSC.hg38", package = "scMORE")
  # Step 2: Initiate GRN object and select candidate regions
  single_cell <- Pando::initiate_grn(single_cell,
                                     peak_assay="peaks",
                                     rna_assay='RNA',
                                     regions=phastConsElements20Mammals.UCSC.hg38,
                                     )

  # Step 3: Scan candidate regions for TF binding motifs
  # Get motif data
  data("motifs", package = "scMORE")
  data("motif2tf", package = "scMORE")
  single_cell <- Pando::find_motifs(
    single_cell,
    pfm = motifs,                          # Position Frequency Matrix (PFM) for motifs
    genome = BSgenome.Hsapiens.UCSC.hg38   # Human genome reference
  )

  # Step 4: Infer the gene regulatory network
  single_cell <- Pando::infer_grn(single_cell,
                                 peak_to_gene_method = peak2gene_method,
                                 method = infer_method)

  # Step 5: Identify regulatory modules
  single_cell <- Pando::find_modules(single_cell,
                                     p_thresh=0.1,
                                     nvar_thresh=2,
                                     min_genes_per_module=1,
                                     rsq_thresh=0.05)

  # Step 6: Extract GRN modules (regulons)
  regulons <- Pando::NetworkModules(single_cell)

  # Step 7: Extract regulatory network data
  data_regulons <- data.frame(
    TF = regulons@meta$tf,                # Transcription factors
    Target = regulons@meta$target,        # Target genes
    Regions = regulons@meta$regions,      # Regulatory regions
    Pval = regulons@meta$pval             # Peak-gene association significance
  )

  # Step 8: Filter regulons with at least `n_genes` target genes
  tf_gene_counts <- data.frame(table(data_regulons$TF))
  valid_tfs <- tf_gene_counts$Var1[tf_gene_counts$Freq >= n_targets]


  # Filter GRN to include only valid TFs
  filtered_regulons <- data_regulons[data_regulons$TF %in% valid_tfs, ]

  # link peaks to target


  # Step 9: Prepare outputs
  grn_outputs <- list(
    grn = filtered_regulons,              # Filtered GRN data
    tf_names = as.vector(valid_tfs)       # List of TFs with valid regulons

  )

  return(grn_outputs)
}
