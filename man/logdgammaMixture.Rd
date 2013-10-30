\name{logdgammaMixture}
\alias{logdgammaMixture}
\title{Compute log density of a gamma mixture}
\usage{
logdgammaMixture(mui, shape, rate, prob)
}
\arguments{
  \item{mui}{density is evaluated at \eqn{\mu_i}}

  \item{shape}{shape parameter \eqn{\alpha} of the gamma
  distribution}

  \item{rate}{rate parameter \eqn{\beta} of the gamma
  distribution}

  \item{prob}{mixture probability \eqn{p}}
}
\value{
  value of log density
}
\description{
  Computes the log density of a mixture of a gamma
  distribution and a point mass at zero
  \eqn{log((1-p)I(\mu_i=0) + p \text{dgamma}(\mu_i, \alpha,
  \beta))}{log((1-p)I(\mu_i=0) + p \text{dgamma}(\mu_i,
  \alpha, \beta))}
}
\keyword{internal}
