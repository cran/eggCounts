\name{fecrtCI}
\alias{fecrtCI}
\title{Compute standard FECRT according to WAAVP guidelines}
\usage{
fecrtCI(epg1, epg2, paired = FALSE, alpha = 0.05, R = 1999)
}
\arguments{
\item{epg1}{faecal egg counts in untreated animals}

\item{epg2}{faecal egg counts in treated animals}

\item{paired}{logical. If TRUE, indicates the samples are paired. Otherwise they are unpaired.}

\item{alpha}{confidence level of the interval}

\item{R}{number of bootstrap replicates}

}
\value{
A list with
\item{estimate}{the estimated percentage reduction in mean epg}
\item{bootCI}{corresponding bootstrap confidence interval}
\item{approxCI}{corresponding approximate confidence interval}
}
\description{
Computes the standard Faecal Egg Count Reduction Test together with approximate
confidence interval according to the WAAVP guidelines
(Coles et al., 1992, 2006). The function also returns bootstrap confidence intervals.
}
\examples{
data(epgs)
fecrtCI(epgs$before, epgs$after, paired=TRUE)
}
\author{
Michaela Paul
}
\references{
Coles GC, Bauer C, Borgsteede FHM, Geerts S, Klei TR, Taylor MA,  Waller, PJ (1992).
World Association for the Advancement of Veterinary Parasitology (WAAVP)
methods for the detection of anthelmintic resistance in nematodes of
veterinary importance, Veterinary Parasitology, 44:35-44.

Coles GC, Jackson F, Pomroy WE, Prichard RK,  von Samson-Himmelstjerna G,
Silvestre A, Taylor MA,  Vercruysse J (2006).
The detection of anthelmintic resistance in nematodes of veterinary	importance,
VeterinaryParasitology, 136:167-185.
}
\keyword{modelling}
