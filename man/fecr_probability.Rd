\name{fecr_probs}
\alias{fecr_probs}
\title{Compute the probability of the reduction parameter relative to the given threshold}
\usage{
fecr_probs(stanFit, threshold = 0.95, lessthan = TRUE, 
           plot = TRUE, xlab, ylab, main, verbose = TRUE, ...)
}
\arguments{
\item{stanFit}{A stanfit object from the output of \code{fecr_stan()}}

\item{threshold}{numeric. The default threshold is 0.95 or 95\%.}

\item{lessthan}{logical. If TRUE, the probability less than the threshold is computed. Otherwise greater or equal to the threshold is computed. Default is TRUE.}

\item{plot}{logical. If TRUE, the posterior density of the reduction is plotted with region less than the threshold shaded.}

\item{xlab, ylab, main}{Arguments for plotting. Only used if \code{showPlot = TRUE}.}

\item{verbose}{logical, If TRUE, a statement with computed probability is printed out.}

\item{...}{Additional plotting arguments.}
}
\value{
Returns a numeric value indicating the probability in percentage.
}
\description{
Computes the probability of the reduction parameter's marginal posterior density relative to a threshold.
}
\author{
Craig Wang
}

\examples{
\dontrun{
## load sample data
data(epgs)

## apply zero-infation model to the data vector
model <- fecr_stan(epgs$before, epgs$after, rawCounts=FALSE, preCF=10, 
                   paired=TRUE, zeroInflation=TRUE)
fecr_probs(model$stan.samples)
}}
\keyword{modelling}

