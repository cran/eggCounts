##################################################################
## Illustrate faecal egg count modelling
##################################################################

# load libraries
library('eggCounts')

##################################################
##  1-sample situation
##################################################
# simulate faecal egg counts
set.seed(123)
counts <- simData1s(n = 20,  # number of samples samples
                    m = 500, # mean 
                    k = 0.5, # overdispersion parameter
                    psi=1,   # prevalence
                    f=50)    # correction factor

# look at simulated counts: matrix with columns 
# (obs = raw counts*correction factor, master= raw counts ,true = true EpG)
head(counts,3)
 
# specify priors and proposal distributions for the 
# one-sample Poisson-Gamma model formulation
priors.mu <- list(hyperpars = c(1, 0.001), proposalDist = "kl")
priors.phi <- list(hyperpars = c(1, 0.7), proposalDist = "unif", v = 0.5)

# use small numbers for testing
nburnin <- 5e3
nsamples <- 1e4
thin <- 1

# run a Poisson-gamma model
resultsP <- fec_mcmc(counts[,"obs"], 
		     rawCounts=FALSE, # r
                     f=50,  # correction factor
                     model="PoGa",   
                     priors.mu=priors.mu, priors.phi=priors.phi, 
		     nburnin=nburnin, nsamples=nsamples, thin=thin)

resultsP

# get samples and covert them to mcmc-objects, so that the functionality 
# of the R-package coda (summary, plot,...) can be used
samplesP <- samples2mcmc(resultsP)
# this is a list with
#  fec  - mean(EpG rate )
#  all1 - parameters of the model except the individual ones
#  mui  - individual epg rates
#  y    - individual true numbers of epg

# for instance
plot(samplesP$all1)
geweke.diag(samplesP$all1)


##################################################
##  2-sample situation
##################################################
# load epgs data before and after treatment
data(epgs)
# plot FECs
epgsL <- reshape(epgs, direction="long",varying=list(c("before","after")))
epgsL$time <- factor(epgsL$time, levels=1:2, labels=c("untreated","after treatment"))

xyplot(before/50 ~ time, group=id, data=epgsL, type=c("p","l"), col=1,
      xlab="", ylab="Faecal egg counts")

# specify priors and proposal distributions
priors.mu <- list(hyperpars = c(1, 0.001))
priors.phi <- list(hyperpars = c(1, 0.7))
priors.delta <- list(hyperpars = c(1,1))

nburnin <- 5e3
nsamples <- 1e4
thin <- 1

## run a paired Poisson-gamma model
result2 <- fecr_mcmc(epgs$before, epgs$after, f.pre=10, f.post=10, model="paired",
           priors.mu=priors.mu, priors.phi=priors.phi, priors.delta=priors.delta,
           nburnin=nburnin, nsamples=nsamples, thin=thin)

result2

# get samples 
samples2 <- samples2mcmc(result2)
# this is a list with
#  fecr   - mean(EpG rate in treated animals), 
#           mean(EpG rate in untreated animals),
#           fecr=reduction in means)
#  all1   - parameters of the model except the individual means
#  muiControl - individual epg rate of control animals
#  muiTreated - individual epg rate of treated animals

# plot samples of model parameters
plot(samples2$all1)
# plot reduction in mean
densityplot(samples2$fecr)

