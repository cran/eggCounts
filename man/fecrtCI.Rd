\name{fecrtCI}
\alias{fecrtCI}
\alias{fecrtCI_Levecke}
\title{Compute standard FECRT according to WAAVP guidelines}
\usage{
fecrtCI(epg1, epg2, paired = FALSE, alpha = 0.05, R = 1999)
fecrtCI_Levecke(epg1, epg2, paired = FALSE, alpha = 0.05, R = 1999)
}
\arguments{
\item{epg1}{numeric vector. Faecal egg counts in untreated animals.}

\item{epg2}{numeric vector. Faecal egg counts in treated animals.}

\item{paired}{logical. If TRUE, indicates samples are paired. Otherwise samples are unpaired.}

\item{alpha}{numeric. Confidence level of the interval.}

\item{R}{numeric. Number of bootstrap replicates.}

}
\value{
A list with
\item{estimate}{estimated percentage reduction in mean epg}
\item{bootCI}{bootstrap confidence interval}
\item{approxCI}{approximate confidence interval}
}
\description{
\code{fecrtCI()} computes the standard Faecal Egg Count Reduction Test
together with approximate confidence interval according to the WAAVP
guidelines (Coles et al., 1992, 2006). 

\code{fecrtCI_Levecke()} computes according to the paper Levecke (2015). 

Both functions also returns bootstrap confidence intervals.
}
\examples{
data(epgs)
fecrtCI(epgs$before, epgs$after, paired = TRUE)
}
\author{
Michaela Paul, Zhixuan Li
}
\references{
Coles GC, Bauer C, Borgsteede FHM, Geerts S, Klei TR, Taylor MA and Waller, PJ (1992).
World Association for the Advancement of Veterinary Parasitology (WAAVP) methods for the detection of anthelmintic resistance in nematodes of veterinary importance, Veterinary Parasitology, 44:35-44.

Coles GC, Jackson F, Pomroy WE, Prichard RK, von Samson-Himmelstjerna G, Silvestre A, Taylor MA and Vercruysse J (2006). 
The detection of anthelmintic resistance in nematodes of veterinary
importance, Veterinary Parasitology, 136:167-185.

Levecke B, et al. (2015).
Chapter Five - Mathematical Inference on Helminth Egg Counts in Stool and Its Applications in Mass Drug Administration Programmes to Control Soil-Transmitted Helminthiasis in Public Health,
in Advances in Parasitology,  Anderson RM, Basanez GM, eds.
87:193-247. DOI:10.1016/bs.apar.2015.01.001

}
\keyword{nonparametric}
