\name{getPrior_delta}
\alias{getPrior_delta}
\title{Get prior parameters from Beta distribution}
\usage{
getPrior_delta(lower, upper, p = 0.7, mode, conc, plot = TRUE)
}
\arguments{
\item{lower,upper,p}{numeric. The lower and upper threshold where \code{p} probability density mass lies in a beta distribution.}

\item{mode,conc}{numeric. The mode and concentration parameters of a beta distribution. Higher concentration indicates smaller variance.}

\item{plot}{logical. If TRUE, the prior distribution is plotted after parameters are found.}
}
\value{
Returns Gamma prior parameters for \eqn{\mu}.
}
\description{
Compute the shape parameters from a Beta distribution for \eqn{\mu} based on some prior belief.
}
\author{
Tea Isler \cr
Craig Wang 
}
\details{
\code{multiroot} function from \strong{rootSolve} package is used to compute the parameters.
}
\examples{
getPrior_delta(lower = 0.6, upper = 0.9, p = 0.8)
}

