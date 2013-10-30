\name{update_phi_unif}
\alias{update_phi_unif}
\title{Update the overdispersion parameter \eqn{phi}
using a uniform random walk proposal}
\usage{
update_phi_unif(phi, mu, sum_mui, sumlog_mui, n, a, 
    b, v.phi)
}
\arguments{
  \item{phi}{phi}

  \item{mu}{mu}

  \item{sum_mui}{sum of \eqn{mu_i}}

  \item{sumlog_mui}{sum of log(\eqn{mu_i})}

  \item{a}{shape parameter of the gamma prior for
  \eqn{phi}}

  \item{b}{rate parameter of the gamma prior for \eqn{phi}}
}
\value{
  updated parameter \eqn{phi}
}
\description{
  Update the overdispersion parameter \eqn{phi} using a
  uniform random walk proposal
}
\keyword{internal}
