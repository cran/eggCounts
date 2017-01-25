\name{fecr_probability}
\alias{fecr_probability}
\title{Compute the probability of the reduction parameter relative to the given threshold}
\usage{
fecr_probability(stanFit, threshold = 0.95, 
  lessthan = TRUE)
}
\arguments{
\item{stanFit}{A stanfit object from the output of fecr_stan()}

\item{threshold}{numeric. The default threshold is 95\% or 0.95.}

\item{lessthan}{logical. If true, the probability less than the threshold is computed. Otherwise greater or equal to the threshold is computed. Default is TRUE.}
}
\value{
Returns the probability in percentage.
}
\description{
Computes the probability of the reduction parameter's marginal posterior density is less than (<), greater or equal to (>=) the threshold.
}
\author{
Craig Wang \email{craig.wang@uzh.ch} 
}

\examples{
\dontrun{
## load the sample data as a vector
data(epgs)

## apply zero-infation model to the data vector
model<-fecr_stan(epgs[,1],epgs[,2],rawCounts=FALSE,preCF=10,paired=TRUE,zeroInflation=TRUE)
fecr_probability(model$stan.samples)
}}


