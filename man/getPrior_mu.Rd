\name{getPrior_mu}
\alias{getPrior_mu}
\title{Get prior parameters from Gamma distribution}
\usage{
getPrior_mu(x, px, y, py, s1 = 1, s2 = 0.001, plot = TRUE)
}
\arguments{
\item{x,px,y,py}{numeric. Threshold of some prior belief about true epg. There is \code{px} probability that the true epg is below \code{x}, and there is \code{py} probability that the true epg is below \code{y}}

\item{s1,s2}{numeric. Starting values.}

\item{plot}{logical. If TRUE, the prior distribution is plotted after parameters are found.}
}
\value{
Returns Gamma prior parameters for \eqn{\mu}.
}
\description{
Compute the shape and rate parameters from a Gamma distribution for \eqn{\mu} based on some prior belief about its cumulative distribution function.
}
\author{
Tea Isler \cr
Craig Wang 
}
\details{
\code{multiroot} function from \strong{rootSolve} package is used to compute the parameters.
}
\examples{
getPrior_mu(x = 200, px = 0.3, y = 500, py = 0.8)
}

