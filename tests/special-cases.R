library(eggCounts)


set.seed(123)

horse2 <- c( 0,0,0,0,0,0,0,0,0,0,0,0,0)
horse1 <- c(50,0,0,0,0,0,0,0,0,0,0,0,0)

p   <- fecr_mcmc(horse1,horse2,nburnin=10,nsamples=250, model="paired",    maxiter.pilot=5, priors.phi=list(v=10), verbose=F)
up  <- fecr_mcmc(horse1,horse2,nburnin=10,nsamples=250, model="unpaired",  maxiter.pilot=5, priors.phi=list(v=10), verbose=F)
u   <- fecr_mcmc(horse1,horse2,nburnin=10,nsamples=250, model="ZIPoGa_u",  maxiter.pilot=5, priors.phi=list(v=10), verbose=F)
upd <- fecr_mcmc(horse1,horse2,nburnin=10,nsamples=250, model="ZIPoGa_upd",maxiter.pilot=5, priors.phi=list(v=10), verbose=F)

print(p,   digits=2, quantiles=c(.1,.25,.5,.75,.9))
print(up,  digits=2, quantiles=c(.1,.25,.5,.75,.9))
print(u,   digits=2, quantiles=c(.1,.25,.5,.75,.9))
print(upd, digits=2, quantiles=c(.1,.25,.5,.75,.9))



p   <- fecr_mcmc(horse2,horse1,nburnin=10,nsamples=250,model="paired",    maxiter.pilot=5, priors.phi=list(v=10), verbose=F)
up  <- fecr_mcmc(horse2,horse1,nburnin=10,nsamples=250,model="unpaired",  maxiter.pilot=5, priors.phi=list(v=10), verbose=F)
u   <- fecr_mcmc(horse2,horse1,nburnin=10,nsamples=250,model="ZIPoGa_u",  maxiter.pilot=5, priors.phi=list(v=10), verbose=F)
upd <- fecr_mcmc(horse2,horse1,nburnin=10,nsamples=250,model="ZIPoGa_upd",maxiter.pilot=5, priors.phi=list(v=10), verbose=F)

print(p,   digits=2, quantiles=c(.1,.25,.5,.75,.9))
print(up,  digits=2, quantiles=c(.1,.25,.5,.75,.9))
print(u,   digits=2, quantiles=c(.1,.25,.5,.75,.9))
print(upd, digits=2, quantiles=c(.1,.25,.5,.75,.9))


