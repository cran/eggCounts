\name{update_psi_beta}
\alias{update_psi_beta}
\title{Update the zero-inflation parameter \eqn{psi}
using a beta approximation to the FC as proposal}
\usage{
update_psi_beta(psi, nNZ_b, nZ_b, nNZ_a, delta, n, 
    a, b)
}
\arguments{
  \item{psi}{psi}

  \item{nNZ_b}{nNZ.b}

  \item{nZ_b}{nZ.b}

  \item{nNZ_a}{nNZ.a}

  \item{delta}{delta}

  \item{n}{n}

  \item{a}{shape1 parameter of the beta prior for
  \eqn{psi}}

  \item{b}{shape2 parameter of the beta prior for
  \eqn{psi}}
}
\value{
  updated parameter \eqn{phi}
}
\description{
  Update the zero-inflation parameter \eqn{psi} using a
  beta approximation to the FC as proposal
}
\keyword{internal}
