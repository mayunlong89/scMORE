
get_gmt <- function(geneset) {
  # Convert each gene set to GMT format as a string
  gmt_data <- sapply(names(geneset), function(i) {
    paste(c(i, "temp", geneset[[i]]), collapse = "\t")
  })

  # Return as a character vector
  return(gmt_data)
}
