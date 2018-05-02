\name{fec_stan}
\alias{fec_stan}
\title{Modelling of faecal egg count data (one-sample case)}
\usage{
fec_stan(fec, rawCounts = FALSE, CF = 50, zeroInflation = TRUE, 
  muPrior, kappaPrior, phiPrior, 
  nsamples = 2000, nburnin = 1000, thinning = 1, nchain = 2, 
  ncore = 1, adaptDelta = 0.95, saveAll = FALSE, verbose = FALSE)
}
\arguments{
\item{fec}{vector of faecal egg counts}

\item{rawCounts}{logical. If true, \code{preFEC} and \code{postFEC} correspond to raw counts
(as counted on equipment). Otherwise they correspond to calculated epgs (raw counts times correction factor).
Defaults to \code{FALSE}.}

\item{CF}{a positive integer or a vector of positive integers. Correction factor(s)}

\item{zeroInflation}{logical. If true, uses the model with zero-inflation. Otherwise uses the model without zero-inflation}

\item{muPrior}{a list with hyper-prior information for the group mean epg parameter \eqn{\mu}. The default prior is \code{list(priorDist = "gamma",hyperpars=c(1,0.001))}, i.e. a gamma distribution with shape 1 and rate 0.001, its 90\% probability mass lies between 51 and 2996}

\item{kappaPrior}{a list with hyper-prior information for the group dispersion parameter \eqn{\kappa}. The default prior is \code{list(priorDist = "gamma",hyperpars=c(1,0.7))}, i.e. a gamma distribution with shape 1 and rate 0.7, its 90\% probability mass lies between 0.1 and 4.3 with a median of 1}

\item{phiPrior}{a list with hyper-prior information for zero-inflation parameter. The default prior is \code{list(priorDist = "beta",hyperpars=c(1,1))}, i.e. a uniform prior between 0 and 1}

\item{nsamples}{a positive integer specifying the number of samples for each chain (including burn-in samples)}

\item{nburnin}{a positive integer specifying the number of burn-in samples}

\item{thinning}{a positive integer specifying the thinning parameter, the period for saving samples}

\item{nchain}{a positive integer specifying the number of chains}

\item{ncore}{a positive integer specifying the number of cores to use when executing the chains in parallel}

\item{adaptDelta}{the target acceptance rate, a numeric value between 0 and 1}

\item{saveAll}{logical. If TRUE, posterior samples for all parameters are saved in the \code{stanfit} object. If FALSE, only samples for \eqn{\mu}, \eqn{\kappa} and \eqn{\phi} are saved. Default to FALSE.}

\item{verbose}{logical. If true, prints progress and debugging information}

}
\value{
Prints out summary of \code{meanEPG} as the posterior mean epg. The posterior summary contains the mean, standard deviation (sd), 2.5\%,  50\% and 97.5\% percentiles, the 95\% highest posterior density interval (HPDLow95 and HPDHigh95) and the posterior mode. NOTE: we recommend to use the 95\% HPD interval and the mode for further statistical analysis.

The returned value is a list that consists of:
\item{stan.samples}{An object of S4 class \code{stanfit} representing the fitted results. For more information, please see the \href{https://cran.r-project.org/package=rstan}{\code{stanfit-class}} in \code{rstan} reference manual.}

\item{posterior.summary}{A data frame that is the same as the printed posterior summary.}
}
\description{
Models faecal egg counts data in a one-sample case with (zero-inflated) Poisson-gamma model formulation using Stan modelling language. It is computationally several-fold faster compared to conventional MCMC techniques. For the installation instruction of Stan, please read  \href{https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started}{https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started}.
}
\details{
The first time each non-default model is applied, it can take up to 20 seconds for Stan to compile the model. Currently the function only support prior distributions with two parameters. For a complete list of supported priors and their parameterization, please consult the list of distributions in Stan \href{http://mc-stan.org/documentation/}{http://mc-stan.org/documentation/}.

The default number of samples per chain is 2000, with 1000 burn-in samples. Normally this is sufficient in Stan. If the chains do not converge, one should tune the MCMC parameters until convergence is reached to ensure reliable results.
}
\author{
Craig Wang 
}
\seealso{
\code{\link{simData1s}} for simulating faecal egg count data with one sample 
}
\examples{
\dontrun{
## load the sample data
data(epgs)

## apply zero-infation model
model <- fec_stan(epgs$before, rawCounts=FALSE, CF=50)
samples <- stan2mcmc(model$stan.samples)
}}
\keyword{modelling}