\name{MH_RW_unif}
\alias{MH_RW_unif}
\title{MH step with uniform (truncated at 0) proposal around the current value}
\usage{
MH_RW_unif(old, logUpost, v)
}
\arguments{
  \item{old}{current state of the chain}

  \item{logUpost}{function evaluating the log unnormalized
  marginal posterior}

  \item{v}{tuning parameter, detemines the width of the
  uniform proposal}
}
\value{
  new state of the chain with a boolean indicating
  acceptance as attribute
}
\description{
  MH step with uniform (truncated at 0) proposal around the
  current value
}
\keyword{internal}
