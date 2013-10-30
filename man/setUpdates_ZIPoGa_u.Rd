\name{setUpdates_ZIPoGa_u}
\alias{setUpdates_ZIPoGa_u}
\title{Specify parameter updates for ZIPoGa_u (2-sample) model}
\usage{
setUpdates_ZIPoGa_u(priors, fec.pre, fec.post, f.pre, 
    f.post)
}
\arguments{
  \item{priors}{list with prior specifications}

  \item{fec.pre}{vector with faecal egg counts (before
  treatment)}

  \item{fec.post}{vector with faecal egg counts (after
  treatment)}

  \item{f.pre}{correction factors before treatment}

  \item{f.post}{correction factors after treatment}
}
\value{
  a named list with update functions
}
\description{
  Specify parameter updates for ZIPoGa_u (2-sample) model
}
\keyword{internal}
