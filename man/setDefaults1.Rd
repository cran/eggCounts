\name{setDefaults1}
\alias{setDefaults1}
\title{Set default values for the one-sample (ZI)PoGa model formulation}
\usage{
setDefaults1(priors.mu, priors.phi, priors.psi)
}
\arguments{
  \item{priors.mu}{named list with hyperprior
  specifications for \eqn{\mu}{mu}, containing elements
  \code{hyperpars} and \code{proposalDist}}

  \item{priors.phi}{named list with hyperprior
  specifications for the overdispersion parameter
  \eqn{\phi}{phi}, containing elements \code{hyperpars},
  \code{proposalDist} (unif or gamma) and in the case of a
  uniform proposal also a tuning parameter \code{v}}

  \item{priors.psi}{named list with hyperprior
  specifications for the prevalence \eqn{\psi}{psi},
  containing elements \code{hyperpars}}
}
\value{
  A named list with prior specifications for \eqn{\mu}{mu},
  \eqn{\phi}{phi} and \eqn{\psi}{psi}.  Default prior
  distributions are: \item{mu}{gamma[1, 0.001]}
  \item{phi}{gamma[1, 0.1]} \item{psi}{beta[1, 1]} Default
  proposal distributions for \eqn{\mu}{mu} and
  \eqn{\phi}{phi} are: \item{mu}{approximating inverse
  gamma, gamma or lognormal distribution}
  \item{phi}{unif[max(current.value - v,0), current.value +
  v], with \eqn{v=0.5}}
}
\description{
  Set default values for the one-sample (ZI)PoGa model
  formulation
}
\keyword{internal}
