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
A MCMC object with a list of relevant parameters depending on the chosen model.
}

\description{
Converts a large \code{stanfit} object into a MCMC object for easier analysis, it extracts the relevant MCMC samples of the model from stanFit including the baseline mean epg, dispersion, pre- and post-treatment zero-inflation parameters and the calculated reduction.
}
\details{
The output can be analyzed as a typical MCMC object with the functions from the \code{coda} package. NOTE: The resulting MCMC object does not contain warm-up samples and is already thinned.  
}
\examples{
\dontrun{
data(epgs)

## apply zero-infation model for the paired design to the data vector
model <- fecr_stan(epgs[,1], epgs[,2], rawCounts=FALSE, preCF=10, paired=TRUE, zeroInflation=TRUE)
samples <- stan2mcmc(model$stan.samples)
summary(samples)
}}
\author{
Craig Wang \email{craig.wang@uzh.ch} 
}

