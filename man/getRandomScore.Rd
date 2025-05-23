% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getRandomScore.R
\name{getRandomScore}
\alias{getRandomScore}
\title{Generate Random Regulon Scores}
\usage{
getRandomScore(
  tf_list,
  target_scores_background,
  real_specificity,
  real_importance,
  len_of_regulon,
  theta = 0.5,
  alpha = 1,
  top_n = 5
)
}
\arguments{
\item{tf_list}{A vector of transcription factors (TFs) to randomly select from.}

\item{target_scores_background}{A data frame containing background scores for target genes.
The data frame must have columns \code{genes} and scores for filtering and sampling.}

\item{real_specificity}{A numeric vector of real specificity scores for sampling.}

\item{real_importance}{A numeric vector of real importance scores for sampling.}

\item{len_of_regulon}{An integer specifying the total number of genes (including the TF) in the regulon.
Must be greater than 1.}

\item{theta}{A numeric value representing the weight for combining TF specificity with gene scores.
Default is 0.5.}

\item{alpha}{A numeric value representing the weight for the penalty term in the final score.
Default is 1.}

\item{top_n}{An integer specifying the number of top genes (ranked by importance) to include in TF calculations.
Default is 5.}
}
\value{
A list containing the following elements:
\item{SpecificityScore}{The calculated specificity score (M1) of the regulon.}
\item{ImportanceWeightScore}{The calculated importance weight score (M2) of the regulon.}
\item{RegulonScore}{The final calculated regulon score, which incorporates a penalty term.}
}
\description{
This function generates random regulon scores for a given set of transcription factors (TFs),
target gene background scores, and real specificity/importance scores. It computes three scores:
SpecificityScore (M1, CTS), ImportanceWeightScore (M2,GRS), and the final RegulonScore (TRS), which includes a penalty term.
}
