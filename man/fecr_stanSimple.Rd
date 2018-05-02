\name{fecr_stanSimple}
\alias{fecr_stanSimple}
\title{Modelling the reduction of faecal egg count data using a simple Bayesian model}
\usage{
fecr_stanSimple(preFEC, postFEC, rawCounts = FALSE, preCF = 50, postCF = preCF, 
  muPrior, deltaPrior, nsamples = 2000, nburnin = 1000, thinning = 1, nchain = 2, 
  ncore = 1, adaptDelta = 0.95, saveAll = FALSE, verbose = FALSE)
}
\arguments{
\item{preFEC}{vector of pre-treatment faecal egg counts}

\item{postFEC}{vector of post-treatment faecal egg counts}

\item{rawCounts}{logical. If TRUE, \code{preFEC} and \code{postFEC} correspond to raw counts
(as counted on equipment). Otherwise they correspond to calculated epgs (raw counts times correction factor).
Defaults to \code{FALSE}.}

\item{preCF}{a positive integer or a vector of positive integers. Pre-treatment correction factor(s)}

\item{postCF}{a positive integer or a vector of positive integers. Post-treatment correction factor(s)}

\item{muPrior}{a list with hyper-prior information for the group mean epg parameter \eqn{\mu}. The default prior is \code{list(priorDist = "gamma",hyperpars=c(1,0.001))}, i.e. a gamma distribution with shape 1 and rate 0.001, its 90\% probability mass lies between 51 and 2996}

\item{deltaPrior}{a list with hyper-prior information for the reduction parameter \eqn{\delta}. The default prior is \code{list(priorDist = "beta",hyperpars=c(1,1))}, i.e. a uniform prior between 0 and 1}

\item{nsamples}{a positive integer specifying the number of samples for each chain (including burn-in samples)}

\item{nburnin}{a positive integer specifying the number of burn-in samples}

\item{thinning}{a positive integer specifying the thinning parameter, the period for saving samples}

\item{nchain}{a positive integer specifying the number of chains}

\item{ncore}{a positive integer specifying the number of cores to use when executing the chains in parallel}

\item{adaptDelta}{the target acceptance rate, a numeric value between 0 and 1}

\item{saveAll}{logical. If TRUE, posterior samples for all parameters are saved in the \code{stanfit} object. If FALSE, only samples for \eqn{\delta}, \eqn{\mu}, \eqn{\kappa} and \eqn{\phi} are saved. Default to FALSE.}

\item{verbose}{logical. If TRUE, prints progress and debugging information}

}
\value{
Prints out the posterior summary of \code{FECR} as the reduction, \code{meanEPG.untreated} as the mean pre-treatment epg, and \code{meanEPG.treated} as the mean after-treatment epg. The posterior summary contains the mean, standard deviation (sd), 2.5\%,  50\% and 97.5\% percentiles, the 95\% highest posterior density interval (HPDLow95 and HPDHigh95) and the posterior mode. 

NOTE: we recommend to use the 95\% equal-tailed credible interval and the median as summary statistics of reduction for the individual efficacy model. For all other models, we recommend to use the 95\% HPD interval and the mode. 

The returned value is a list that consists of:
\item{stan.samples}{An object of S4 class \code{stanfit} representing the fitted results. For more information, please see the \href{https://cran.r-project.org/package=rstan}{\code{stanfit-class}} in \code{rstan} reference manual.}

\item{posterior.summary}{A data frame that is the same as the printed posterior summary.}
}
\description{
Models the reduction in faecal egg counts data with simple Bayesian model formulation using Stan modelling language. The model is for paired design only, and it assumes analytical sensitivity adjusted Poisson distribution for the observed egg counts. 

It is computationally several-fold faster compared to conventional MCMC techniques. For the installation instruction of Stan, please read: \href{https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started}{Stan Installation.}
}
\details{
The first time each model with non-default priors is applied, it can take up to 20 seconds for Stan to compile the model. Currently the function only support prior distributions with two parameters. For a complete list of supported priors and their parameterization, please consult the list of distributions in \href{http://mc-stan.org/documentation/}{Stan}. 

The default number of samples per chain is 2000, with 1000 burn-in samples. Normally this is sufficient in Stan. If the chains do not converge, one should tune the MCMC parameters until convergence is reached to ensure reliable results.
}

\author{
Tea Isler
Craig Wang 
}
\seealso{
\code{\link{simData2s}} for simulating faecal egg counts data with two samples 
}
\examples{
\dontrun{
## load sample data
data(epgs)

## apply paired model with individual efficacy
model <- fecr_stanSimple(epgs$before, epgs$after, rawCounts=FALSE, preCF=10)
samples <- stan2mcmc(model$stan.samples)

}}
\keyword{modelling}

