\name{fecr_mcmc}
\alias{fecr_mcmc}
\title{Modelling of Faecal Egg Count data (two-sample case)}
\usage{
fecr_mcmc(fec.pre, fec.post, rawCounts = FALSE, f.pre = 50, 
    f.post = f.pre, model = "paired", priors.mu = list(hyperpars = c(1, 
        0.001), proposalDist = "kl"), priors.phi = list(hyperpars = c(1, 
        0.1), proposalDist = "unif", v = 0.5), priors.delta = list(priorDist = c("gamma", 
        "beta")[1], hyperpars = c(1, 1), proposalDist = NULL), 
    priors.psiB = list(hyperpars = c(1, 1)), priors.psiA = list(hyperpars = c(1, 
        1)), priors.deltaPsi = list(hyperpars = c(1, 1), proposalDist = "beta"), 
    maxiter.pilot = 15, nburnin = 5000, nsamples = 10000, thin = 1, 
    initials = NULL, verbose = TRUE, ...)
}
\arguments{
  \item{fec.pre}{vector with faecal egg counts before
  treatment}

  \item{fec.post}{vector with faecal egg counts after
  treatment}

  \item{rawCounts}{logical indicating whether
  \code{fec.pre} and \code{fec.post} corresponds to raw
  counts (as counted on the McMaster slide), or to
  calculated EpGs (raw counts times correction factor).
  Defaults to \code{FALSE}.}

  \item{f.pre}{correction factor(s) before treatment}

  \item{f.post}{correction factor(s) after treatment}

  \item{model}{string with model formulation ("paired" or
  "unpaired")}

  \item{priors.mu}{list with hyper-prior/proposal
  information for \eqn{\mu}{mu}}

  \item{priors.phi}{list with hyper-prior/proposal
  information for \eqn{\phi}{phi}}

  \item{priors.delta}{list with hyper-prior/proposal
  information for \eqn{\delta}{delta}}

  \item{priors.psiA}{list with hyper-prior/proposal
  information for \eqn{\psi_A}{psiA}}

  \item{priors.psiB}{list with hyper-prior/proposal
  information for \eqn{\psi_B}{psiB}}

  \item{priors.deltaPsi}{list with hyper-prior/proposal
  information for \eqn{\delta_\psi}{deltaPsi}}

  \item{maxiter.pilot}{maximal number of tries to determine
  a good tuning value for the proposal distribution for
  \eqn{\phi}{phi}}

  \item{nburnin}{number of burn-in iterations}

  \item{nsamples}{number of desired samples}

  \item{thin}{thinning parameter}

  \item{initials}{named list with starting values for the
  parameters \code{mu}, \code{phi}, \code{muiPre},
  \code{muiPost}, \code{yPre}, \code{yPost}, \code{delta},
  \code{psiB} and \code{psiA}}

  \item{verbose}{print progress information}

  \item{\dots}{extra arguments (not used atm)}
}
\value{
  a named list with \item{samples}{ list with samples and
  acceptance rates} \item{priors}{ list with prior
  specifications} \item{v.phi}{ tuning parameter for
  \eqn{phi}} \item{initials}{ list with initial values}
  \item{model}{ name of the specified model}
  \item{nburnin}{ number of burnin iterations}
  \item{nsamples}{ number of returned samples} \item{thin}{
  used thinning factor}
}
\description{
  Modelling of Faecal Egg Count data in a two-sample case
  using a paired or unpaired (ZI) Poisson gamma model
  formulation
}
\seealso{
  \code{demo("fecm", package = "eggCounts")}
}
