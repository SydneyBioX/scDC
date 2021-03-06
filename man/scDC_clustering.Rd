% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scDC_clustering.R
\name{scDC_clustering}
\alias{scDC_clustering}
\title{scDC_clustering}
\usage{
scDC_clustering(
  exprsMat = NULL,
  cellTypes = NULL,
  subject = NULL,
  calCI = TRUE,
  calCI_method = c("BCa", "multinom", "percentile"),
  nboot = NULL,
  conf_level = 0.95,
  ncores = 1,
  verbose = TRUE
)
}
\arguments{
\item{exprsMat}{logcounts expression matrix with each row represents gene, and each column represents cell}

\item{cellTypes}{A vector indicates the cell type info of the data}

\item{subject}{A vector indicates the subject info of the data}

\item{calCI}{A logical input for whether calculating the confidence interval for proportion}

\item{calCI_method}{A string indicates the method that is used to calculate confidence interval. Options include \code{BCa}, \code{percentile}, and \code{multinom}.}

\item{nboot}{Number of bootstrap. If \code{calCI = TRUE}, \code{nboot = 10000} by default. Otherwise, \code{nboot = 500}.}

\item{conf_level}{confidence level, with default 0.95}

\item{ncores}{Number of cores that are used.}

\item{verbose}{A logical input for whether print the progress.}
}
\value{
Returns a data frame.
}
\description{
Single-cell Differential Composition Analysis with performing clustering
}
\examples{
## Loading example data
library(scDC)
data("sim")

cellTypes = sim$sim_cellTypes
subject = sim$sim_subject
\dontrun{
res_noCALCI = scDC_clustering(cellTypes, subject, calCI = FALSE)
res_BCa = scDC_clustering(cellTypes, subject, calCI = TRUE, calCI_method = "BCa")
res_percentile = scDC_clustering(cellTypes, subject, calCI = TRUE, calCI_method = "percentile")
res_multinom = scDC_clustering(cellTypes, subject, calCI = TRUE, calCI_method = "multinom")

}

}
\author{
Yingxin Lin
}
