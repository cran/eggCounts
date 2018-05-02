\name{stan2mcmc}
\alias{stan2mcmc}
\title{Convert a Stanfit object to MCMC object}
\usage{
stan2mcmc(stanFit)
}
\arguments{
\item{stanFit}{A stanfit object from the output of either \code{fecr_stan()} or \code{fec_stan()}}
}
\value{
A MCMC object with a list of relevant parameters.
}

\description{
Converts a \code{stanfit} object into a MCMC object for easier analysis.
}
\details{
The output can be analyzed as a typical MCMC object with the functions from the \code{coda} package. NOTE: The resulting MCMC object does not contain warm-up samples and is already thinned.  
}
\examples{
\dontrun{
data(epgs)

## apply zero-infation model for the paired design 
model <- fecr_stan(epgs$before, epgs$after, rawCounts = FALSE, 
                   preCF = 10, paired = TRUE, zeroInflation = TRUE)
samples <- stan2mcmc(model$stan.samples)
summary(samples)
}}
\author{
Craig Wang
}
\keyword{modelling}
