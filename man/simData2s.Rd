\name{simData2s}
\alias{simData2s}
\title{Simulate faecal egg count data (2-sample situation)}
\usage{
simData2s(n = 10, preMean = 500, delta = 0.1, kappa = 0.5, 
  phiPre = 1, phiPost = phiPre, f = 50, paired = TRUE, 
  rounding = TRUE)
}
\arguments{
\item{n}{sample size (number of faecal samples collected pre- and post-treatment)}

\item{preMean}{true number of eggs per gram (epg) (i.e. worm burden) before treatment}

\item{delta}{proportion of epg left after treatment, between 0 and 1. 1 - \eqn{\delta} is reduction in mean after treatment, \code{delta} = 0.1 indicates a 90\% reduction}

\item{kappa}{overdispersion parameter, \eqn{\kappa \to \infty}{\kappa -> \infty} corresponds to Poisson}

\item{phiPre}{pre-treatment prevalence (i.e. proportion of infected animals), between 0 and 1}

\item{phiPost}{post-treatment prevalence, between 0 and 1}

\item{f}{correction factor of the egg counting technique, either an integer or a vector of integers with length \code{n}}

\item{paired}{logical. If true, paired samples are simulated. Otherwise unpaired samples are simulated.}

\item{rounding}{logical. If true, the Poisson mean for the raw counts is rounded. The rounding applies since the mean epg is frequently reported as an integer value. For more information, please see Details.}

}
\value{
A matrix with six columns, namely the observed epg (\code{obs}),
    number of eggs counted on microscope slide (\code{master}) and true egg counts (\code{true}) for both pre- and post- treatment.
}
\description{
Generates two samples of (zero-inflated) egg count data
}
\details{
The simulation process does not exactly match the proposed models in [ref:paper], however the simulated data resembles the data observed in real world. 

In the simulation of raw (\code{master}) counts, it follows a Poisson distribution with some mean. The mean is frequently rounded down if it has a very low value and \code{rounding = TRUE},there expects to be a up to 3-10\% positive bias in the mean reduction when \eqn{\mu} < 150 and \eqn{\delta} < 0.1. Set \code{rounding = FALSE} if one does not wish to have any bias.
}
\author{
Craig Wang \email{craig.wang@uzh.ch} \cr
Michaela Paul 
}
\examples{
fec <- simData2s(n=10, preMean=500, delta=0.1, kappa=0.5)

## show the positive bias when the true reduction should be 95\%
set.seed(1)
fec <- simData2s(n=1e5, preMean=150, delta=0.05, kappa=0.5)
1-mean(fec[,5])/mean(fec[,2])
}
\seealso{
\code{\link{fecr_stan}} for analyzing faecal egg count data with two samples 
}

