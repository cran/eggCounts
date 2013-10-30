\name{setUpdates_ZIPoGa2}
\alias{setUpdates_ZIPoGa2}
\title{Specify parameter updates for PoGa2 (1-sample) model}
\usage{
setUpdates_ZIPoGa2(priors, fec, f)
}
\arguments{
  \item{priors}{list with priors for \eqn{\mu}{mu},
  \eqn{\phi}{phi} and \eqn{\psi}{psi}}

  \item{fec}{vector with faecal egg counts}

  \item{f}{correction factor (e.g. 50 for the McMaster
  technique)}
}
\value{
  a named list with update functions
}
\description{
  Specify parameter updates for PoGa2 (1-sample) model
}
\keyword{internal}
