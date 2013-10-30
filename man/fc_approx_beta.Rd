\name{fc_approx_beta}
\alias{fc_approx_beta}
\title{Approximate a full conditional with given mode and curvature by a beta distribution}
\usage{
fc_approx_beta(mode, curvature, dens01)
}
\arguments{
  \item{mode}{mode of the full conditional}

  \item{curvature}{curvature at the mode}

  \item{dens01}{density at the boundary}
}
\value{
  shape1 and shape2 parameter of the approximating beta
  distribution
}
\description{
  Approximate a full conditional with given mode and
  curvature by a beta distribution
}
\keyword{internal}
