% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convertRegulon.R
\name{convertRegulon}
\alias{convertRegulon}
\title{Convert scMORE TRS Results Format for Calculating Energy Score}
\usage{
convertRegulon(regulon2disease_results, targetCelltype = 1)
}
\arguments{
\item{regulon2disease_results}{TRS results from regulon2disease() function.}

\item{targetCelltype}{Can be specified as either the index of the target cell type column (e.g., 1),
or the name of the target cell type (e.g., "Monocytes").}
}
\description{
Transform scMORE TRS results for the getDistance() function by manually choosing a target cell type.
The TRS values of other cell types are averaged as the non-target cell type.
}
