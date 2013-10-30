\name{update_h_invGamma}
\alias{update_h_invGamma}
\title{Update the parameter \eqn{theta} (with FC of form h) using an inverse gamma approximation to the FC as proposal}
\usage{
update_h_invGamma(theta, argA, argC, a, b)
}
\arguments{
  \item{argA}{argument a (in parts) of the function h (e.g.
  n*phi for theta=mu)}

  \item{argC}{argument c of the function h (e.g.
  sum(mui)*phi for theta=mu)}

  \item{a}{shape parameter of the gamma prior for
  \eqn{theta}}

  \item{b}{rate parameter of the gamma prior for
  \eqn{theta}}
}
\value{
  updated parameter \eqn{theta}
}
\description{
  Update the parameter \eqn{theta} (with FC of form h)
  using an inverse gamma approximation to the FC as
  proposal
}
\keyword{internal}
