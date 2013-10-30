\name{update_theta_gammaMix}
\alias{update_theta_gammaMix}
\title{Update a parameter vector theta \eqn{theta_1,\ldots,theta_n} based on a gamma mixture distribution}
\usage{
update_theta_gammaMix(mui, y, mu, phi, psi, n, whichZero_fec, 
    nNonZero_fec, p = 1)
}
\arguments{
  \item{mui}{mui}

  \item{y}{y}

  \item{mu}{mu}

  \item{phi}{phi}

  \item{psi}{mixture probability}

  \item{n}{length of theta}

  \item{nNonZero_fec}{number of non-zero observed counts}

  \item{whichZero_fec}{vector of logicals indicating which
  observed counts are zero}

  \item{p}{p, defaults to 1}
}
\value{
  updated parameter \eqn{theta}
}
\description{
  Update a parameter vector theta
  \eqn{theta_1,\ldots,theta_n} based on a gamma mixture
  distribution
}
\keyword{internal}
