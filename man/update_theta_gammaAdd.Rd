\name{update_theta_gammaAdd}
\alias{update_theta_gammaAdd}
\title{Update a parameter vector theta \eqn{theta_1,\ldots,theta_n} based on a gamma distribution
using a gamma proposal with slightly shifted shape parameter to avoid to small theta}
\usage{
update_theta_gammaAdd(theta, s, r, n)
}
\arguments{
  \item{theta}{theta}

  \item{s}{shape}

  \item{r}{rate}

  \item{n}{length of theta}
}
\value{
  updated parameter \eqn{theta}
}
\description{
  Update a parameter vector theta
  \eqn{theta_1,\ldots,theta_n} based on a gamma
  distribution using a gamma proposal with slightly shifted
  shape parameter to avoid to small theta
}
\keyword{internal}
