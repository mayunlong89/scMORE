% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getRandomScore.R
\name{getRandomScore}
\alias{getRandomScore}
\title{Calculate the random scores for background genes}
\usage{
getRandomScore(target_scores_sub)
}
\arguments{
\item{target_scores_sub}{A data frame containing specificity scores and MAGMA z-scores for each gene.
Must contain the following columns:
- scores: Specificity scores for each gene.
- magma_zscore: MAGMA z-scores for each gene.}
}
\value{
A data frame with the original gene names and generated random specificity and MAGMA z-scores.
}
\description{
Generate random scores based on the fitted distributions of specificity scores and MAGMA z-scores.
}
