#' @title Link SNPs to genomic peaks
#'
#' @param snp_info Data frame containing SNP information (required columns: CHR, POS, LP, SNP)
#' @param peak2gene_strength List containing GRN output, with a 'Regions' column for peaks
#' @param buffer Numeric value specifying the flanking sizes (default = 0)
#'        - extend peak range (upstream and downstream + 50bp)
#'        - buffer = 50bp, represent add 50 bp to both upstream and downstream of peak
#' @return Data frame mapping SNPs to genomic peaks with associated log-transformed p-values
#' @importFrom GenomicRanges GRanges resize start end seqnames mcols
#' @importFrom IRanges IRanges width findOverlaps
#' @export
#'
snp2peak <- function(snp_info, peak2gene_strength, buffer = 50) {
  # Step 1: Format SNP bed file
  # Ensure required columns are present and chromosome format is consistent
  snp_info <- checkSNPs(snp_info, add_chr_prefix = TRUE)

  # Prepare SNP bed format
  snp_bed <- snp_info %>%
    dplyr::mutate(end = POS + 1) %>%
    dplyr::rename(chr = CHR, start = POS, log_P = LP, snp_id = SNP) %>%
    dplyr::select(chr, start, end, log_P, snp_id)

  # Ensure chromosome names have "chr" prefix
  snp_bed <- getSNPchr(snp_bed)

  # Create SNP GenomicRanges object
  snp_ranges <- GenomicRanges::GRanges(
    seqnames = snp_bed$chr,
    ranges = IRanges::IRanges(start = snp_bed$start, end = snp_bed$end),
    snp_id = snp_bed$snp_id,
    log_P = snp_bed$log_P
  )

  # Step 2: Validate and format peak regions
  if (!"Regions" %in% colnames(peak2gene_strength)) {
    stop("peak2gene_strength must contain a 'Regions' column")
  }

  peak2gene_strength$peak_ids <- peak2gene_strength$Regions
  # Process Regions column
  peak_bed <- peak2gene_strength %>%
    dplyr::mutate(Regions = as.character(Regions)) %>%
    tidyr::separate(Regions, into = c("chr", "start", "end"), sep = "-", fill = "right") %>%
    dplyr::mutate(
      start = as.numeric(start),
      end = as.numeric(end)
    ) %>%
    dplyr::filter(!is.na(start) & !is.na(end)) %>%
    dplyr::select(chr, start, end, TF, Target, Strength, peak_ids)

  # Create Peak GenomicRanges object
  peak_ranges <- GenomicRanges::GRanges(
    seqnames = peak_bed$chr,
    ranges = IRanges::IRanges(start = peak_bed$start, end = peak_bed$end),
    TF = peak_bed$TF,
    Target = peak_bed$Target,
    peak_ids = peak_bed$peak_ids,
    Strength = peak_bed$Strength
  )

  # Extend peak ranges
  if (buffer > 0) {
    extended_peak_ranges <- GenomicRanges::resize(
      peak_ranges,
      width = IRanges::width(peak_ranges) + buffer * 2,
      fix = "center"
    )
    GenomicRanges::start(extended_peak_ranges) <- pmax(GenomicRanges::start(extended_peak_ranges), 1) # Ensure ranges start >= 1
  } else {
    extended_peak_ranges <- peak_ranges
  }
  GenomicRanges::mcols(extended_peak_ranges) <- GenomicRanges::mcols(peak_ranges)

  # Step 3: Find overlaps
  overlaps <- IRanges::findOverlaps(snp_ranges, extended_peak_ranges)

  # Step 4: Extract overlap information
  overlap_results <- data.frame(
    snp_id = GenomicRanges::mcols(snp_ranges)$snp_id[queryHits(overlaps)],
    peak_chr = GenomicRanges::seqnames(extended_peak_ranges)[subjectHits(overlaps)],
    peak_start = GenomicRanges::start(extended_peak_ranges)[subjectHits(overlaps)],
    peak_end = GenomicRanges::end(extended_peak_ranges)[subjectHits(overlaps)],
    peak_ids = GenomicRanges::mcols(extended_peak_ranges)$peak_ids[subjectHits(overlaps)],
    TF = GenomicRanges::mcols(extended_peak_ranges)$TF[subjectHits(overlaps)],
    Target = GenomicRanges::mcols(extended_peak_ranges)$Target[subjectHits(overlaps)],
    Strength = GenomicRanges::mcols(extended_peak_ranges)$Strength[subjectHits(overlaps)],
    Significance = GenomicRanges::mcols(snp_ranges)$log_P[queryHits(overlaps)]
  )

  # Step 5: Ensure all peaks are included
  all_peaks <- data.frame(
    peak_chr = GenomicRanges::seqnames(extended_peak_ranges),
    peak_start = GenomicRanges::start(extended_peak_ranges),
    peak_end = GenomicRanges::end(extended_peak_ranges),
    peak_ids = GenomicRanges::mcols(extended_peak_ranges)$peak_ids,
    TF = GenomicRanges::mcols(extended_peak_ranges)$TF,
    Target = GenomicRanges::mcols(extended_peak_ranges)$Target,
    Strength = GenomicRanges::mcols(extended_peak_ranges)$Strength
  )

  # Merge to include non-overlapping peaks
  final_results <- dplyr::full_join(
    all_peaks, overlap_results,
    by = c("peak_chr", "peak_start", "peak_end", "peak_ids", "TF", "Target", "Strength")
  )

  return(final_results)
}
