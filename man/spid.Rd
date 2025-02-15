% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SL_wilcoxon.R
\name{spid}
\alias{spid}
\title{perform Wilcoxon for potential geneA - geneB synthetic lethal pair
check if geneA altered cell lines are more dependent on geneB}
\usage{
spid(
  sl_pair,
  sens_matrix1,
  alt_matrix1,
  plot = F,
  experiment_name = "",
  min_alt_no = 10
)
}
\arguments{
\item{sl_pair}{c(geneA_ENSEMBL, geneB_ENSEMBL) from a list of potential SL gene pairs}

\item{min_alt_no}{minimum no of cell lines which have geneA alteration, default 10}

\item{sens_matrix}{geneB dependency matrix in CCLE cell lines for a given knock-down experiment}

\item{alt_matrix}{geneA alteration matrix for CCLE cell lines}
}
\value{
p-value of the Wilcoxon test
}
\description{
perform Wilcoxon for potential geneA - geneB synthetic lethal pair
check if geneA altered cell lines are more dependent on geneB
}
