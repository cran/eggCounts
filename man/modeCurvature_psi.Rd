\name{modeCurvature_psi}
\alias{modeCurvature_psi}
\title{Compute mode and curvature at the mode of a full conditional of form
\eqn{h(x) = (1-d x)^{a} (1-x)^{b} x^{c}}{h(x) = (1-d*x)^a (1-x)^b x^c}}
\usage{
modeCurvature_psi(a, b, c, d)
}
\arguments{
  \item{a}{a}

  \item{b}{b}

  \item{c}{c}

  \item{d}{d}
}
\value{
  list with mode, curvature at the mode, and the density
  const*h(mode) should the mode be at the boundary
}
\description{
  Compute mode and curvature at the mode of a full
  conditional of form \eqn{h(x) = (1-d x)^{a} (1-x)^{b}
  x^{c}}{h(x) = (1-d*x)^a (1-x)^b x^c}
}
\keyword{internal}
