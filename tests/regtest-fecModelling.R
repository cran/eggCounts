library(eggCounts)

# simulate egg counts
set.seed(123)
counts <- simData1s(n = 10, m = 500, k = 0.5)

# test different dilution factors
set.seed(1)
f1 <- fec_mcmc(counts[,"obs"], rawCounts=FALSE, f=50, model=("PoGa"), nburnin=10,nsamples=1e3,thin=1, maxiter.pilot=1, verbose=F)
set.seed(1)
f2 <- fec_mcmc(counts[,"obs"], rawCounts=FALSE, f=rep(50,10), model=("PoGa"), nburnin=10,nsamples=1e3,thin=1, maxiter.pilot=1, verbose=F)
set.seed(1)
f3 <- fec_mcmc(counts[,"obs"]/50, rawCounts=TRUE, f=rep(50,10), model=("PoGa"), nburnin=10,nsamples=1e3,thin=1, maxiter.pilot=1, verbose=F)

stopifnot(all.equal(f1, f2))
stopifnot(all.equal(f1, f3))
