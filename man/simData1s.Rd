\name{simData1s}
\alias{simData1s}
\title{Simulate faecal egg count data (1-sample situation)}
\usage{
simData1s(n = 10, m = 500, k = 0.5, psi = 1, f = 50)
}
\arguments{
  \item{n}{sample size (number of faeces collected)}

  \item{m}{true number of eggs per gram (EPG) (i.e. worm
  burden)}

  \item{k}{overdispersion factor, \eqn{k \to \infty}{k ->
  Inf} corresponds to Poisson}

  \item{psi}{prevalence (i.e. proportion of infected
  animals)}

  \item{f}{correction factor of the egg counting technique}
}
\value{
  a data-frame with the observed EPG (\code{obs}), number
  of eggs counted on microscope slide (\code{master}) and
  true egg counts (\code{true}).
}
\description{
  Simulates (zero-inflated) egg count data
}
\examples{
fec <- simData1s(n = 10, m = 500, k = 0.5, psi = 0.7)
}