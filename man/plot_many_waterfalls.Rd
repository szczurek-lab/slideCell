% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/waterfalls_Wilcoxon.R
\name{plot_many_waterfalls}
\alias{plot_many_waterfalls}
\title{creates 4 waterfall plots
where geneA altered cell lines (for Wilcoxon test subpopulation) are colored
and geneB depenency is plotted}
\usage{
plot_many_waterfalls(
  sl_pair,
  sens_matrix = list(SLiCell::sens_Drive_RSA, SLiCell::sens_Drive_Ataris,
    SLiCell::sens_Achilles_AC, SLiCell::sens_Achilles_Demeter),
  data_set_names = c("Drive_RSA", "Drive_Ataris", "Achilles_AC", "Achilles_Demeter"),
  alt_matrix1 = SLiCell::gene_cl_mut_matrix
)
}
\arguments{
\item{sl_pair}{list of geneA geneB pairs for verifing}

\item{sens_matrix}{list of geneB dependency matrices in geneA altered cell lines, default - all experimental data preloaded in the package}

\item{data_set_names}{names for the dependency data sets - will be visible in results}

\item{alt_matrix1}{geneA alteration matrix for cell lines, default - mutation_matrix preloaded in the package}
}
\value{
4 waterfall plots
}
\description{
creates 4 waterfall plots
where geneA altered cell lines (for Wilcoxon test subpopulation) are colored
and geneB depenency is plotted
}
