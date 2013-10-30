\name{fc_approx_gamma}
\alias{fc_approx_gamma}
\title{Approximate a full conditional with given mode and curvature by a gamma distribution}
\usage{
fc_approx_gamma(mode, curvature)
}
\arguments{
  \item{mode}{mode of the full conditional}

  \item{curvature}{curvature at the mode}
}
\value{
  shape and rate parameter of the approximating gamma
  distribution
}
\description{
  Approximate a full conditional with given mode and
  curvature by a gamma distribution
}
\keyword{internal}
