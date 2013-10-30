\name{update_theta_beta}
\alias{update_theta_beta}
\title{Update a parameter \eqn{theta} based on a beta distribution}
\usage{
update_theta_beta(x, n, a, b)
}
\arguments{
  \item{x}{number of events}

  \item{n}{sample size}

  \item{a}{shape1 parameter of the beta prior for
  \eqn{theta}}

  \item{b}{shape2 parameter of the beta prior for
  \eqn{theta}}
}
\value{
  updated parameter \eqn{theta}
}
\description{
  Update a parameter \eqn{theta} based on a beta
  distribution
}
\keyword{internal}
