\name{setInitials2}
\alias{setInitials2}
\title{Set initial values for all parameters of the two-sample (ZI)PoGa model}
\usage{
setInitials2(initials, fec.pre, fec.post, f.pre, f.post)
}
\arguments{
  \item{initials}{named list with initial values for
  parameters \code{mu} \code{phi} \code{muiPre}
  \code{muiPost} \code{yPre} \code{yPost} \code{delta}
  \code{psiB} \code{psiA} and \code{deltaPsi}}

  \item{fec.pre}{faecal egg count data before treatment}

  \item{fec.post}{faecal egg count data after treatment}

  \item{f.pre}{correction factors before treatment}

  \item{f.post}{correction factors after treatment}
}
\value{
  A named list with initial values
}
\description{
  Set initial values for all parameters of the two-sample
  (ZI)PoGa model
}
\keyword{internal}
