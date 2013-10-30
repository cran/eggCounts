
##' Simulate faecal egg count data (1-sample situation)
##'
##' Simulates (zero-inflated) egg count data 
##'
##' @param n sample size (number of faeces collected)
##' @param m true number of eggs per gram (EPG) (i.e. worm burden)
##' @param k overdispersion factor, \eqn{k \to \infty}{k -> Inf} corresponds to Poisson
##' @param psi prevalence (i.e. proportion of infected animals)
##' @param f correction factor of the egg counting technique
##' @return a data-frame with the observed EPG (\code{obs}), 
##'     number of eggs counted on microscope slide (\code{master}) and 
##'     true egg counts (\code{true}).
##' @export

##' @examples
##' fec <- simData1s(n=10, m=500, k=0.5, psi=0.7)
simData1s <- function(n=10, m=500, k=0.5, psi=1, f=50){

   # sample true egg counts
   mean.true <- rgamma(n,shape=k, rate=k/m)
   infected <- sample(c(FALSE,TRUE),size=n, replace=TRUE,prob=c(1-psi,psi))
   mean.true[!infected] <- 0
   epg.true <- rpois(n,mean.true)
   # take a subsample and count eggs y, such that epg = f*y
   fec <- rpois(n,lambda=round(epg.true/f))

   data <- cbind(obs=fec*f, master=fec, true=epg.true)
   return(data)
}


##' Simulate faecal egg count data (2-sample situation)
##'
##' Generates two samples of (zero-inflated) egg count data 
##'
##' @param n sample size (number of faeces collected pre and post treatment)
##' @param m.pre true number of eggs per gram (EPG) (i.e. worm burden) before treatment
##' @param delta reduction in mean after treatment
##' @param k overdispersion factor, \eqn{k \to \infty}{k -> Inf} corresponds to Poisson
##' @param psi prevalence (i.e. proportion of infected animals) 
##' @param psiA prevalence after treatment 
##' @param f correction factor of the egg counting technique
##' @param paired logical indicating a paired or unpaired situation
##' @return a data-frame with the observed EPGs (\code{obs.prePost}), 
##'     number of eggs counted on microscope slide (\code{master}) and true egg counts (\code{true}).
##' @export

##' @examples
##' fec <- simData2s(n=10, m.pre=500, delta=0.8, k=0.5, psi=0.7)
simData2s <- function(n=10, m.pre=500, delta=0.1, k=0.5,psi=1,psiA=psi, f=50, paired=TRUE){

  if(paired){
	# sample true egg counts before treatment
	mean.true <- rgamma(n,shape=k, rate=k/m.pre)
	infected <- sample(c(FALSE,TRUE),size=n, replace=TRUE,prob=c(1-psi,psi))
	mean.true[!infected] <- 0
	epg.true <- rpois(n,mean.true)
	# take a subsample and count eggs y, such that epg = f*y
	fec <- rpois(n,lambda=round(epg.true/f))

	# now add sample egg counts after treatment
	infectedA <- infected
	infectedA[infected] <- sample(c(FALSE,TRUE),size=sum(infected), replace=TRUE,prob=c(1-psiA,psiA))

	mean.true.post <- mean.true*delta
	mean.true.post[!infectedA] <- 0
	epg.true.post <- rpois(n,mean.true.post)
	# take a subsample and count eggs y, such that epg = f*y
	fec.post <- rpois(n,lambda=round(epg.true.post/f))

	data <- cbind(obs.pre =fec*f, master.pre=fec, true.pre=epg.true, obs.post=fec.post*f, master.post=fec.post, true.post=epg.true.post)

  } else {
	data.pre <- simData1s(n=n, m=m.pre, k=k, psi=psi)
	data.post <- simData1s(n=n, m=m.pre*delta, k=k, psi=psiA)

	data <- cbind(data.pre,data.post)
    colnames(data) <- paste(colnames(data),rep(c("pre","post"),each=3),sep=".")
	  
  }

   return(data)

}

