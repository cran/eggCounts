\name{update_h_gamma}
\alias{update_h_gamma}
\title{Update the parameter \eqn{theta}
(with FC of form \eqn{h(x) = x^{-a} \exp(-bx -c/x)}{h(x) = x^(-a) exp(-bx -c/x)})
using a gamma approximation to the FC as proposal}
\usage{
update_h_gamma(theta, argA, argC, a, b)
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
  Update the parameter \eqn{theta} (with FC of form
  \eqn{h(x) = x^{-a} \exp(-bx -c/x)}{h(x) = x^(-a) exp(-bx
  -c/x)}) using a gamma approximation to the FC as proposal
}
\keyword{internal}
