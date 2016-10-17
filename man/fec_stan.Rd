\name{fec_stan}
\alias{fec_stan}
\title{Modelling of faecal egg count data (one-sample case)}
\usage{
fec_stan(fec, rawCounts = FALSE, CF = 50, zeroInflation = TRUE, 
  muPrior, kappaPrior, phiPrior, nsamples = 12000, 
  nburnin = 2000, thinning = 1, nchain = 1, 
  ncore = 1, adaptdelta = 0.9, verbose = FALSE)
}
\arguments{
\item{fec}{vector of faecal egg counts}

\item{rawCounts}{logical. If true, \code{preFEC} and \code{postFEC} correspond to raw counts
(as counted on the McMaster slide). Otherwise they correspond to calculated epgs (raw counts times correction factor).
Defaults to \code{FALSE}.}

\item{CF}{correction factor or vector of correction factors}

\item{zeroInflation}{logical. If true, uses the model with zero-inflation. Otherwise uses the model without zero-inflation}

\item{muPrior}{a list with hyper-prior information for the baseline mean parameter \eqn{\mu}. The default prior is \code{list(priorDist = "gamma",hyperpars=c(1,0.001))}, i.e. a gamma distribution with shape 1 and rate 0.001, its 90\% probability mass lies between 51 and 2996}

\item{kappaPrior}{a list with hyper-prior information for the dispersion parameter \eqn{\kappa}. The default prior is \code{list(priorDist = "gamma",hyperpars=c(1,0.7))}, i.e. a gamma distribution with shape 1 and rate 0.7, its 90\% probability mass lies between 0.1 and 4.3 with a median of 1}

\item{phiPrior}{a list with hyper-prior information for zero-inflation parameter. The default prior is \code{list(priorDist = "beta",hyperpars=c(1,1))}}

\item{nsamples}{a positive integer specifying how many iterations for each chain (including burn-in samples)}

\item{nburnin}{number of burn-in samples}

\item{thinning}{thinning parameter, a positive integer specifying the period for saving samples}

\item{nchain}{a positive integer specifying the number of chains}

\item{ncore}{number of cores to use when executing the chains in parallel}

\item{adaptdelta}{the target acceptance rate, a value between 0 and 1}

\item{verbose}{logical. If true, prints progress and debugging information}

}
\value{
An object of S4 class \code{stanfit} representing the fitted results. For more information, please see the \href{https://cran.r-project.org/package=rstan}{\code{stanfit-class}} in \code{rstan} reference manual.

Prints out summary of \code{meanEPG} as the posterior mean egg count. The posterior summary contains the mean, standard deviation (sd), 2.5\%, 25\%, 50\%, 75\% and 97.5\% percentiles, the 95\% highest posterior density interval (HPDLow95 and HPDHigh95) and the posterior mode. NOTE: we recommend to use the 95\% HPD interval and the mode for further statistical analysis.
}
\description{
Models faecal egg counts data in a one-sample case with (zero-inflated) Poisson-gamma model formulation using the Stan modelling language. It is computationally several-fold faster compare to conventional MCMC techniques. For the installation instruction of Stan, please read: \href{https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Windows}{Stan Installation.}
}
\details{
The first time each non-default model is applied, it can take up to 20 seconds for stan to compile the model. Currently the function only support prior distributions with two parameters. For a complete list of supported priors and their parameterization, please consult the list of distributions in \href{http://mc-stan.org/documentation/}{Stan}.

Sometimes the function outputs informational message from Stan regarding the Metropolis proposal rejections, this is due to the sampler hit the boundary of the parameter space. For some variables, the boundary point is not supported in the distribution. This is not a concern if there are only a few such warnings. 
}
\author{
Craig Wang \email{craig.wang@uzh.ch} 
}
\seealso{
\code{\link{simData1s}} for simulating faecal egg count data with one sample 
}
\examples{
\dontrun{
## load the sample data as a vector
data(echinococcus)
fec<-echinococcus[[1]]

## apply zero-infation model to the data vector
model<-fec_stan(fec,rawCounts=FALSE,CF=50,zeroInflation=FALSE)
samples<-stan2mcmc(model)

## a demonstration
demo("fecm_stan", package = "eggCounts") 
}}
