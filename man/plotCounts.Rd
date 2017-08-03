\name{plotCounts}
\alias{plotCounts}
\title{Plot eggCounts data}
\usage{
plotCounts(data, paired=TRUE, points=TRUE, points.method="jitter", 
           xlabel="", ylabel="Faecal egg counts", ...)
}
\arguments{
\item{data}{a data frame with two columns, the first column is before treatment counts, the second column is after treatment counts.}

\item{paired}{logical. If true, uses the plot for the paired design. Otherwise uses the plot for the unpaired design}

\item{points}{logical. If true, add individual points for unpaired plot. Omitted if paired is TRUE.}

\item{points.method}{a quoted keyword to be used to separate coincident points if points is TRUE. The default method "overplot" causes such points to be overplotted, but it is also possible to specify "jitter" to jitter the points, or "stack" have coincident points stacked.}

\item{xlabel}{label of x-axis.}

\item{ylabel}{label of y-axis.}

\item{...}{Additional argument for function \code{xyplot} if paired is TRUE, for function \code{boxplot} otherwise.}
}
\value{
A plot is returned based on the arguments.
}
\description{
Plot egg count data to reflect changes between before and after treatment.
}
\author{
Craig Wang \email{craig.wang@uzh.ch} \cr
}
\details{
For paired data, a xyplot is used. For unpiared data, a grouped boxplot is used.
}
\examples{
data(epgs)
plotCounts(epgs)
}

