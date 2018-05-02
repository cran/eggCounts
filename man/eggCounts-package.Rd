\docType{package}
\name{eggCounts-package}
\alias{eggCounts}
\alias{eggCounts-package}
\title{Hierarchical modelling of faecal egg counts}
\description{
This package implements Bayesian hierarchical models for the analysis of faecal egg count data. Bayesian inference is done via efficient MCMC sampling using Stan. Additional models are available externally for handling FECs with potential outliers or bi-modality. The models are in \strong{eggCountsExtra} package hosted on Github. 
}
\details{
\tabular{ll}{
Package: \tab eggCounts\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2018-04-30\cr
License: \tab GPL (>= 3) \cr
LazyLoad: \tab yes\cr
}}

\author{
Craig Wang \email{craig.wang@uzh.ch} \cr
Michaela Paul 
}

\examples{
\dontrun{

## Citations
citation('eggCounts')

## History of changes
file.show(system.file("NEWS", package = "eggCounts"))

## Demonstration
demo("fecm_stan", package = "eggCounts") 

## Install eggCountsExtra
devtools::install_github("CraigWangUZH/eggCountsExtra")
}}

\keyword{package}

