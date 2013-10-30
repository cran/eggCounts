\name{modeCurvature_phi}
\alias{modeCurvature_phi}
\title{Numerically compute mode and curvature at the mode of a full conditional for \eqn{\phi}{phi} of form
-(n*phi+a-1)*log(phi) + n*lgamma(phi) + n*phi*log(mu)- (phi-1)*logprodmu + phi*(summu.mu+b)}
\usage{
modeCurvature_phi(n, mu, logprodmu, summu.mu, a, b)
}
\arguments{
  \item{n}{n}

  \item{mu}{mu}

  \item{logprodmu}{logprodmu}

  \item{summu.mu}{summu.mu}

  \item{a}{a}

  \item{b}{b}
}
\value{
  list with mode and curvature at the mode
}
\description{
  Numerically compute mode and curvature at the mode of a
  full conditional for \eqn{\phi}{phi} of form
  -(n*phi+a-1)*log(phi) + n*lgamma(phi) + n*phi*log(mu)-
  (phi-1)*logprodmu + phi*(summu.mu+b)
}
\keyword{internal}
