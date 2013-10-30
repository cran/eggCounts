\name{simData2s}
\alias{simData2s}
\title{Simulate faecal egg count data (2-sample situation)}
\usage{
simData2s(n = 10, m.pre = 500, delta = 0.1, k = 0.5, 
    psi = 1, psiA = psi, f = 50, paired = TRUE)
}
\arguments{
  \item{n}{sample size (number of faeces collected pre and
  post treatment)}

  \item{m.pre}{true number of eggs per gram (EPG) (i.e.
  worm burden) before treatment}

  \item{delta}{reduction in mean after treatment}

  \item{k}{overdispersion factor, \eqn{k \to \infty}{k ->
  Inf} corresponds to Poisson}

  \item{psi}{prevalence (i.e. proportion of infected
  animals)}

  \item{psiA}{prevalence after treatment}

  \item{f}{correction factor of the egg counting technique}

  \item{paired}{logical indicating a paired or unpaired
  situation}
}
\value{
  a data-frame with the observed EPGs (\code{obs.prePost}),
  number of eggs counted on microscope slide
  (\code{master}) and true egg counts (\code{true}).
}
\description{
  Generates two samples of (zero-inflated) egg count data
}
\examples{
fec <- simData2s(n = 10, m.pre = 500, delta = 0.8, 
    k = 0.5, psi = 0.7)
}
