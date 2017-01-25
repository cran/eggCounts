\name{simData1s}
\alias{simData1s}
\title{Simulate faecal egg count data (1-sample situation)}
\usage{
simData1s(n = 10, mean = 500, kappa = 0.5, phi = 1, f = 50, rounding = TRUE)
}
\arguments{
\item{n}{sample size (number of faeces collected)}

\item{mean}{true number of eggs per gram (epg)}

\item{kappa}{overdispersion parameter, \eqn{\kappa \to \infty}{\kappa -> \infty} corresponds to Poisson}

\item{phi}{prevalence i.e. proportion of infected animals, between 0 and 1}

\item{f}{correction factor of the egg counting technique, either an integer or a vector of integers with length \code{n}}

\item{rounding}{logical. If true, the Poisson mean for the raw counts is rounded. The rounding applies since the mean epg is frequently reported as an integer value. For more information, please see Details.}
}
\value{
A matrix with three columns, namely the observed epg (\code{obs}),
    number of eggs counted on the McMaster slide (\code{master}) and
    true egg counts (\code{true}).
}
\description{
Simulates (zero-inflated) egg count data
}
\author{
Craig Wang \email{craig.wang@uzh.ch} \cr
Michaela Paul 
}
\details{
The simulation process does not exactly match the proposed models in [ref:paper], however the simulated data resembles the data observed in real world. 

In the simulation of raw (\code{master}) counts, it follows a Poisson distribution with some mean. The mean is frequently rounded down if it has a very low value and \code{rounding = TRUE}, hence there expects to be a negative bias overall when \eqn{\mu} < 150. Set \code{rounding = FALSE} if one does not wish to have any bias in the simulated counts.
}
\examples{
fec <- simData1s(n=10, mean=500, kappa=0.5, phi=0.7)
}
\seealso{
\code{\link{fec_stan}} for analyzing faecal egg count data with one sample 
}

