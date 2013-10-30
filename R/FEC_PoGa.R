#############################################################
##  Faecal Egg Count Modelling (1-sample situation)
#############################################################

##' @include MHsteps.R
##' @include updates.R
##' @include simData.R
{}

##' Set default values for the one-sample (ZI)PoGa model formulation
##'
##' @param priors.mu named list with hyperprior specifications for \eqn{\mu}{mu}, 
##'    containing elements \code{hyperpars} and \code{proposalDist} 
##' @param priors.phi named list with hyperprior specifications for the overdispersion parameter \eqn{\phi}{phi}, 
##'    containing elements \code{hyperpars}, \code{proposalDist} (unif or gamma) and 
##'    in the case of a uniform proposal also a tuning parameter \code{v} 
##' @param priors.psi named list with hyperprior specifications for the prevalence \eqn{\psi}{psi}, 
##'    containing elements \code{hyperpars} 

##' @return A named list with prior specifications for \eqn{\mu}{mu}, \eqn{\phi}{phi}
##'    and \eqn{\psi}{psi}.
##'    Default prior distributions are:
##'    \item{mu}{gamma[1, 0.001]}
##'    \item{phi}{gamma[1, 0.1]}
##'    \item{psi}{beta[1, 1]}
##'   Default proposal distributions for \eqn{\mu}{mu} and \eqn{\phi}{phi} are:
##'    \item{mu}{approximating inverse gamma, gamma or lognormal distribution}
##'    \item{phi}{unif[max(current.value - v,0), current.value + v], with \eqn{v=0.5}}
##'
##' @keywords internal
setDefaults1 <- function(priors.mu, priors.phi, priors.psi){
   if(missing(priors.mu)) priors.mu = list(hyperpars=c(1,0.001), proposalDist="kl")
   if(is.null(priors.mu[["hyperpars", exact = TRUE]])) priors.mu$hyperpars = c(1,0.001)
   if(is.null(priors.mu[["proposalDist", exact = TRUE]])) priors.mu$proposalDist = "kl"
   
   if(missing(priors.phi)) priors.phi = list(hyperpars=c(1,0.1), proposalDist="unif", v=0.5)
   if(is.null(priors.phi[["hyperpars", exact = TRUE]])) priors.phi$hyperpars = c(1,0.1)
   if(is.null(priors.phi[["proposalDist", exact = TRUE]])) priors.phi$proposalDist = "unif"
   if(is.null(priors.phi[["v", exact = TRUE]])) priors.phi$v = 0.5

   if(missing(priors.psi)) priors.psi = list(hyperpars=c(1,1))
   if(is.null(priors.psi[["hyperpars", exact = TRUE]])) priors.psi$hyperpars = c(1,1)

#check that (a,b)-values are positive
   
   return(list(mu=priors.mu,phi=priors.phi, psi = priors.psi))
}

##' Set initial values for all parameters of the one-sample (ZI)PoGa model
##'
##' @param initials named list with initial values for parameters 
##'       \code{mu}  
##'       \code{phi}  
##'       \code{mui}  
##'       \code{y}  
##'       and \code{psi}  
##' @param fec faecal egg count data
##' @param f correction factor

##' @return A named list with initial values

##' @keywords internal
setInitials1 <- function(initials, fec, f){
   # set initials
   y <- fec*f
   mu <-  mean(y)
   mui <- y+rbinom(length(fec),size=1,prob=.5)
   phi <- 0.5 
   psi <- max(0.1, mean(fec==0))

   if(is.null(initials)) initials <- list(mu=mu, phi=phi, y=y, mui=mui, psi=psi)

   if(is.null(initials[["y", exact = TRUE]])) initials$y  <- y
   if(is.null(initials[["mu", exact = TRUE]])) initials$mu  <- mu
   if(is.null(initials[["mui", exact = TRUE]])) initials$mui  <- mui
   if(is.null(initials[["phi", exact = TRUE]])) initials$phi  <- phi
   if(is.null(initials[["psi", exact = TRUE]])) initials$psi  <- psi

   # ToDo: check user-defined initials
   
   return(initials)
}

##' Print information about the progress of the MCMC run
##'
##' @param piter percentage of interations already done
##' @param time.start starting time of the run

##' @return prints progress information about the estimated time left
##' @keywords internal
statusMsg <- function(piter, time.start){
      #Sys.sleep(.2)
      time.current <- Sys.time()
      time.elapsed <- difftime(time.current, time.start, units="secs")
      time.remaining <- (100-piter)/piter*time.elapsed
      units <- ifelse(time.remaining<60,"secs","mins")
      if(piter==0){
	 cat("\nStarted to run the MCMC algorithm at",format(time.start,"%T"),"\b.\n")
      } else if(piter==1){
#      cat(sprintf("Estimated duration: %5.1f %s\n",as.numeric(time.remaining,units=units),units))
      cat(sprintf("Completed: %2.0f%%, time remaining: %5.1f %s",1,as.numeric(time.remaining,units=units),units))
      } else if(piter<100) {
	 cat(sprintf("\rCompleted: %2.0f%%, time remaining: %5.1f %s", piter,as.numeric(time.remaining,units=units),units))
      } else cat("\rCompleted: 100%                   at",format(time.current,"%T"),"\b. \n")
}

##' Modelling of Faecal Egg Count data (one-sample case)
##'
##' Modelling of Faecal Egg Count data in a one-sample case using
##' either a Poisson gamma or a zero-inflated Poisson gamma model
##' formulation.
##'
##' @rdname fec_mcmc
##' @param fec vector with faecal egg counts
##' @param rawCounts logical indicating whether \code{fec} corresponds to raw counts
##'  (as counted on the McMaster slide), or to calculated EpGs (raw counts times correction factor). 
##'  Defaults to \code{FALSE}.
##' @param f correction factor for the McMaster technique (e.g. 50). Either a number
##'  or a vector with different correction factors for each FEC
##' @param model either "PoGa" or "ZIPoGa"
##' @param priors.mu list with hyper-prior/proposal information for \eqn{\mu}{mu}
##' @param priors.phi list with hyper-prior/proposal information for \eqn{\phi}{phi}
##' @param priors.psi list with hyper-prior information for \eqn{\psi}{psi}
##' @param maxiter.pilot maximal number of tries to determine a good tuning value
##' for the proposal distribution for \eqn{\phi}{phi}
##' @param nburnin number of burn-in iterations
##' @param nsamples number of desired samples
##' @param thin thinning parameter
##' @param initials named list with starting values for the parameters 
##'         \code{mu}, \code{phi}, \code{mui}, \code{y},\code{psi}
##' @param verbose print progress information
##' @param \dots extra arguments (not used)

##' @return a named list with 
##' \item{samples}{ list with samples and acceptance rates}
##' \item{priors}{ list with prior specifications}
##' \item{v.phi}{ tuning parameter for \eqn{phi}}
##' \item{initials}{ list with initial values}
##' \item{model}{ name of the specified model}
##' \item{nburnin}{ number of burnin iterations}
##' \item{nsamples}{ number of returned samples}
##' \item{thin}{ used thinning factor}

##' @seealso \code{demo("fecm", package = "eggCounts")}

##' @export
fec_mcmc <- function(
   fec,
   rawCounts = FALSE,
   f = 50, 
   model = c("PoGa","ZIPoGa")[1],
   priors.mu = list(hyperpars= c(1,0.001),
 		    proposalDist = "kl"),
   priors.phi = list(hyperpars= c(1,0.1),
 		     proposalDist = "unif",
		     v = 0.5),
   priors.psi = list(hyperpars= c(1,1)),
   maxiter.pilot = 10,
   nburnin = 1e3,
   nsamples = 1e4,
   thin =1,
   initials = NULL,
   verbose=TRUE,
   ...
   ){

   # check class of fec
   if(inherits(fec,"data.frame")) fec <- as.matrix(fec)
   if(is.vector(fec)) fec <- as.matrix(fec) 
   if(!is.matrix(fec)) stop("Faecal egg counts must be provided as a data.frame, matrix or vector\n")

   if(ncol(fec)!=1) stop("Only 1 column is allowed for the fec data\n")

   # number of faecal samples
   n <- length(fec)

   # check correction factors
   if(length(f)>1 && length(f)!=n) stop("Lengths of the vectors f and fec do not match\n")

   # raw counts or EpGs?
   dilution <- f
   if(rawCounts){
      dilution <- 1
   } 

   # divide data by correction factor
   if(any(fec %% dilution !=0)) stop("Correction factor does not match the given fec\n")
   fec <- fec/dilution

   # check model and set default values
   model <- match.arg(model, c("PoGa","ZIPoGa"))
   priors <- setDefaults1(priors.mu=priors.mu, priors.phi=priors.phi, priors.psi=priors.psi)
 
   # set update functions
   updates <- switch(model, 
                     "PoGa"= setUpdates_PoGa_mu(priors, fec,f),
                     "ZIPoGa"= setUpdates_ZIPoGa_mu(priors, fec, f),
                     "PoGa2" = setUpdates_PoGa2(priors, fec, f),
                     "ZIPoGa2"= setUpdates_ZIPoGa2(priors, fec, f),
               )

   update_y <- updates$y
   update_mui <- updates$mui
   update_mu <- updates$mu
   update_phi <- updates$phi
   update_psi <- updates$psi
   runMCMC <- eval(updates$runMCMC)
   
   # set initials 
   start <- setInitials1(initials, fec, f) 
   v.phi <- priors$phi$v

   ############################################################
   ## determine tuning parameter for phi if required
   ###########################################################
   if(priors$phi$proposalDist %in% c("unif","unif01") & maxiter.pilot>0){
      blocksize <- 500

      phi.AR.ok <- FALSE
      startpilot <- start

      # run pilot chain to determine tuning parameter for phi
      i <- 1
      if(verbose) {
	 cat("\nStarted to run pilot chain at",format(Sys.time(),"%T"),"\b\n")
	 cat("using tuning parameter: v =",v.phi,"\n")
      }
      while(!phi.AR.ok && i <= maxiter.pilot){

	 pilot <- runMCMC(startpilot,nsamples=blocksize,nburnin=100,thin=1, v.phi=v.phi, time.start=NULL, verbose=FALSE)
	 startpilot <- pilot$last
	 AR <- checkAR(pilot$accept["phi"],v.phi, verbose=verbose)
	 v.phi <- AR[1]
	 phi.AR.ok <- AR[2]

	 i <- i+1
      }      
      if(verbose) {
	 cat("Finished to run pilot chain at",format(Sys.time(),"%T"),"\b.\n")
	 cat(sprintf("Selected tuning parameter: v = %.1f (acc rate = %3.1f%%)\n",v.phi,100*pilot$accept["phi"]))
      }
      if(!phi.AR.ok) cat("NOTE: Tuning parameter for phi possibly not ok.\n")
   }
   ########################################################################
   # do updates
   time.start <- Sys.time()
   samples <- runMCMC(start,nsamples=nsamples,nburnin=nburnin,thin=thin, v.phi=v.phi, time.start=time.start, verbose=verbose)

   #return results
   result <- list(samples=samples, priors=priors, v.phi=v.phi,initials=start, model=model, nburnin=nburnin, nsamples=nsamples,thin=thin)
   class(result) <- "fecm"
   return(result)

}

#######################################################
### parameter updates for model PoGa_mu (1-sample)
#######################################################

##' Specify parameter updates for PoGa (1-sample) model
##'
##' @param priors  list with priors for \eqn{\mu}{mu} and \eqn{\phi}{phi}
##' @param fec vector with faecal egg counts
##' @param f correction factor for the McMaster technique (e.g. 50)

##' @return a named list with update functions
##'
##' @keywords internal
setUpdates_PoGa_mu <- function(priors, fec, f){ 

   n <- length(fec)

   #hyperparameters for mean mu
   a.mu <- priors$mu$hyperpars[1]
   b.mu <- priors$mu$hyperpars[2]
   #hyperparameters  for overdispersion parameter phi
   a.phi <- priors$phi$hyperpars[1]
   b.phi <- priors$phi$hyperpars[2]
   v.phi <- priors$phi$v

#check that (a,b)-values are positive
   
   #####################################
   ## set up update functions
   ####################################
   # insert data/hyperparameters as default arguments

   ### latent true egg counts y_i
   update_y <- update_Y
   formals(update_y)$n <- n
   formals(update_y)$propC <- 1-1/f
   formals(update_y)$fec <- fec

   ## mean epg rates
   update_mui <- update_theta_gamma
   #update_mui <- update_theta_gammaAdd
   formals(update_mui)$n <- n

   ### mean mu
   # determine the proposal distribution  
   update_mu <- switch(priors$mu$proposalDist,
			"gamma" = update_h_gamma,
			"invGamma" = update_h_invGamma,
			"logNormal" = update_h_logNormal,
			"kl" = update_h_kl
			)   
   formals(update_mu)$a  <- a.mu
   formals(update_mu)$b <- b.mu

   ### overdispersion parameter phi
   update_phi <- switch(priors$phi$proposalDist,
			"gamma" = update_phi_gamma,
			"unif" = update_phi_unif,
			"unif01" = update_phiTransf_unif
			)      
   formals(update_phi)$a  <- a.phi
   formals(update_phi)$b <- b.phi
   formals(update_phi)$n <- n

   ######
   runMCMC <- bquote(function(start, nsamples,nburnin,thin, v.phi,time.start=NULL, verbose=TRUE, n=.(n)){ 
      # number of iterations
      niter <- nsamples*thin
      # print progress information at the following iteration numbers if time.start is given
      if(!verbose | is.null(time.start)){
	 printMsg <- function(...) {NULL}
      } else printMsg <- statusMsg
      niter.1percent <- c(-nburnin,floor(quantile((-nburnin):niter,prob=1:100/100)))

      ## create empty results 
      y.samples <- mui.samples  <- array(dim=c(nsamples, n))  
      mu.samples <- phi.samples <- numeric(nsamples)
      #accept <- rep(0,n+2)
      #names(accept) <- c("mu","phi", paste("mui",1:n,sep="")) 
      accept <- c("mu"=0,"phi"=0)

      ## set initial Values
      y <- y.samples[1,] <- start$y
      mu <-  mu.samples[1] <- start$mu
      mui <- mui.samples[1,] <- start$mui
      phi <- phi.samples[1] <- start$phi 

      sum_mui <- sum(mui.samples[1,])
      sumlog_mui <- sum(log(mui.samples[1,]))

      # update parameters
      for (i in (-nburnin):niter){
	 # update mean
	 mu <- update_mu(mu, n*phi, sum_mui*phi)   
	 
	 # update latent epg counts
	 y <- update_y(mui)

	 # update true latent means
	 mui <- update_mui(s=y+phi, r=1+phi/mu)
	 #mui <- update_mui(mui,s=y+phi, r=1+phi/mu)
	 sum_mui <- sum(mui)
	 sumlog_mui <- sum(log(mui))

	 # update overdispersion
	 phi <- update_phi(phi, mu, sum_mui, sumlog_mui,v.phi=v.phi)

	 # save samples
	 if(i>0 && i%%thin ==0){
	    idx <- i/thin
	    y.samples[idx,] <- y
	    mui.samples[idx,] <- mui
	    mu.samples[idx] <- mu
	    phi.samples[idx] <- phi
	 }

	 # save acceptance
	 accept <- accept + c(attributes(mu)$accept,attributes(phi)$accept)
	 #accept <- accept + c(attributes(mu)$accept,attributes(phi)$accept,attributes(mui)$accept)

	 # print status information
	 if(length(w <- which(niter.1percent %in% i))) {printMsg(w-1, time.start)}
#cat("\nmu=",mu," phi=",phi, " y=",y, " mui=",mui, "sum_mui=",sum_mui,"\n")
      }
      # save last value of chain
      last <- list(y=y, mu=mu, mui=mui,phi=phi)
      # return results
      res <- list(mu=mu.samples, y=y.samples, mui=mui.samples, phi=phi.samples, accept=accept/(nburnin+1+niter),last=last)
      return(res)
   }
   )
   ####

   result <- list(y=update_y,mui=update_mui, mu=update_mu, phi=update_phi, runMCMC=runMCMC) 
   return(result)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#######################################################
### parameter updates for model PoGa_mu (1-sample, zero-inflated)
#######################################################

##' Specify parameter updates for PoGa (1-sample) model
##'
##' Specify parameter updates for PoGa (1-sample) model
##'
##' @param priors  list with priors for \eqn{\mu}{mu}, \eqn{\phi}{phi} and \eqn{\psi}{psi}
##' @param fec vector with faecal egg counts
##' @param f correction factor (e.g. 50 for the McMaster technique)

##' @return a named list with update functions
##'
##' @keywords internal
setUpdates_ZIPoGa_mu <- function(priors, fec, f){ 

   n <- length(fec)
   # which observed counts are zero
   whichZero_fec <- fec==0
   nNonZero_fec <- sum(fec > 0)
   propC <- 1- 1/f

   #hyperparameters for mean mu
   a.mu <- priors$mu$hyperpars[1]
   b.mu <- priors$mu$hyperpars[2]
   #hyperparameters  for overdispersion parameter phi
   a.phi <- priors$phi$hyperpars[1]
   b.phi <- priors$phi$hyperpars[2]
   v.phi <- priors$phi$v
   # hyperparameters for zero-inflation parameter psi
   a.psi <- priors$psi$hyperpars[1]
   b.psi <- priors$psi$hyperpars[2]

   
   ##########################################
   ## set up update functions
   ##########################################
   # insert data/hyperparameters as default arguments

   #######################################################
   ## latent true counts
   #####################################################
   update_y <- update_Y
   formals(update_y)$n <- n
   formals(update_y)$propC <- 1-1/f
   formals(update_y)$fec <- fec
   
  
   ##################################################
   ## population mean mu
   ####################################################
   update_mu <- switch(priors$mu$proposalDist,
			"gamma" = update_h_gamma,
			"invGamma" = update_h_invGamma,
			"logNormal" = update_h_logNormal,
			"kl" = update_h_kl
			) 
   formals(update_mu)$a <- a.mu
   formals(update_mu)$b <- b.mu

   
   #############################################
   ## overdispersion parameter phi
   ##############################################
   update_phi <- switch(priors$phi$proposalDist,
			"gamma" = update_phi_gamma,
			"unif" = update_phi_unif,
			"unif01" = update_phiTransf_unif
			)      
   formals(update_phi)$a  <- a.phi
   formals(update_phi)$b <- b.phi

   #######################################################
   ## latent means mu_i
   #####################################################
   update_mui <- update_theta_gammaMix
   formals(update_mui)$n <- n
   formals(update_mui)$whichZero_fec <- whichZero_fec
   formals(update_mui)$nNonZero_fec <- nNonZero_fec

   #######################################################
   ## zero-inflation parameter
   #####################################################
   update_psi <- update_theta_beta
   formals(update_psi)$n <- n
   formals(update_psi)$a <- a.psi
   formals(update_psi)$b <- b.psi

   ######
   runMCMC <- bquote(function(start, nsamples,nburnin, thin, v.phi,time.start=NULL, verbose=TRUE, n=.(n), nZero_fec=.(sum(fec==0))){
      # number of iterations
      niter <- nsamples*thin
     # print progress information at the following iteration numbers if time.start is given
      if(!verbose | is.null(time.start)){
	 printMsg <- function(...) {NULL}
      } else printMsg <- statusMsg
      niter.1percent <- c(-nburnin,floor(quantile((-nburnin):niter,prob=1:100/100)))

      ## create empty results 
      y.samples <- mui.samples  <- array(dim=c(nsamples, n))
      mu.samples <- phi.samples <- psi.samples <- numeric(nsamples)
      accept <- numeric(2+nZero_fec)
      names(accept) <- c("mu","phi",paste("mui",1:nZero_fec,sep=".")) 
      ## set initial Values
      y <- y.samples[1,] <- start$y
      mu <-  mu.samples[1] <- start$mu
      mui <- mui.samples[1,] <- start$mui
      phi <- phi.samples[1] <- start$phi 
      psi <- psi.samples[1] <- start$psi

      nNonZero_mui <- sum(mui.samples[1,]>0)
      whichNonZero_mui <- mui.samples[1,] > 0
      sum_mui <- sum(mui.samples[1,])
      sumlog_mui <- sum(log(mui.samples[1,whichNonZero_mui]))

      # update parameters
      for (i in (-nburnin):niter){
	 # update mean
	 mu <- update_mu(mu, nNonZero_mui*phi, phi*sum_mui)   
#cat("\nmu=",mu,"\n")
if(is.na(mu)) break
	 # update latent epg counts
	 y <- update_y(mui)
#cat("y=c(",paste(y,","),"\n")

	 # update true latent means
	 mui <- update_mui(mui,y,mu,phi,psi)
#cat("mui=c(",paste(mui,","),"\n")
	 nNonZero_mui <- sum(mui>0)
	 whichNonZero_mui <- mui > 0
	 sum_mui <- sum(mui)
	 sumlog_mui <- sum(log(mui[whichNonZero_mui]))
   
	 # update overdispersion
	 phi <- update_phi(phi, mu, sum_mui, sumlog_mui,nNonZero_mui, v.phi=v.phi)
#cat("phi=",phi,"\n")

	 # update zero-inflation
	 psi <- update_psi(nNonZero_mui)
#cat("psi=",psi,"\n")

	 # save samples
	 if(i>0 && i%%thin ==0){
	    idx <- i/thin
	    y.samples[idx,] <- y
	    mui.samples[idx,] <- mui
	    mu.samples[idx] <- mu
	    phi.samples[idx] <- phi
	    psi.samples[idx] <- psi
	 }

	 # save acceptance
	 accept <- accept + c(attributes(mu)$accept,attributes(phi)$accept, attributes(mui)$accept)

	 # print status information
	 if(length(w <- which(niter.1percent %in% i))) {printMsg(w-1, time.start)}
#cat("\nmu=",mu," phi=",phi, " psi=",psi, " y=",y, " mui=",mui," nNonZero_mui=",nNonZero_mui, "sum_mui=",sum_mui,"\n")
      }
      # save last value of chain
      last <- list(y=y, mu=mu, mui=mui,phi=phi,psi=psi)
      # return results
      res <- list(mu=mu.samples, y=y.samples, mui=mui.samples, phi=phi.samples,psi=psi.samples, accept=accept/(nburnin+1+niter),last=last)
      return(res)
   }
   )
   ####

   result <- list(y=update_y,mui=update_mui, mu=update_mu, phi=update_phi, psi=update_psi, runMCMC=runMCMC) 
   return(result)
}


##' Compute log density of a gamma mixture
##'
##' Computes the log density of a mixture of a gamma distribution and a point mass at zero
##' \eqn{log((1-p)I(\mu_i=0) + p \text{dgamma}(\mu_i, \alpha, \beta))}{log((1-p)I(\mu_i=0) + p \text{dgamma}(\mu_i, \alpha, \beta))}
##'
##' @param mui density is evaluated at \eqn{\mu_i}
##' @param shape shape parameter \eqn{\alpha} of the gamma distribution
##' @param rate rate parameter \eqn{\beta} of the gamma distribution
##' @param prob mixture probability \eqn{p}

##' @return value of log density
##'
##' @keywords internal
logdgammaMixture <- function(mui, shape, rate, prob){
   ifelse(mui==0,
                log(1-prob),
                log(prob)+dgamma(mui, shape=shape,rate=rate,log=T)
         )   
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# try another formulation, where the stages 
#  obs|true ~ Bin(true,p)  and  true|mui ~Pois(mui)  are combined to
#  obs|mui ~ Pois(p*mui)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##' Specify parameter updates for PoGa2 (1-sample) model
##'
##' @param priors  list with priors for \eqn{\mu}{mu} and \eqn{\phi}{phi}
##' @param fec vector with faecal egg counts
##' @param f correction factor for the McMaster technique (e.g. 50)

##' @return a named list with update functions
##'
##' @keywords internal
setUpdates_PoGa2 <- function(priors, fec, f){ 

   n <- length(fec)

   #hyperparameters for mean mu
   a.mu <- priors$mu$hyperpars[1]
   b.mu <- priors$mu$hyperpars[2]
   #hyperparameters  for overdispersion parameter phi
   a.phi <- priors$phi$hyperpars[1]
   b.phi <- priors$phi$hyperpars[2]
   v.phi <- priors$phi$v

#check that (a,b)-values are positive
   
   #####################################
   ## set up update functions
   ####################################
   # insert data/hyperparameters as default arguments

   ## mean epg rates
   update_mui <- update_theta_gamma
   formals(update_mui)$n <- n

   ### mean mu
   # determine the proposal distribution  
   update_mu <- switch(priors$mu$proposalDist,
			"gamma" = update_h_gamma,
			"invGamma" = update_h_invGamma,
			"logNormal" = update_h_logNormal,
			"kl" = update_h_kl
			)   
   formals(update_mu)$a  <- a.mu
   formals(update_mu)$b <- b.mu

   ### overdispersion parameter phi
   update_phi <- switch(priors$phi$proposalDist,
			"gamma" = update_phi_gamma,
			"unif" = update_phi_unif
			)      
   formals(update_phi)$a  <- a.phi
   formals(update_phi)$b <- b.phi
   formals(update_phi)$n <- n

   ######
   runMCMC <- bquote(function(start, nsamples,nburnin,thin, v.phi,time.start=NULL, verbose=TRUE, n=.(n), prop=.(1/f), fec=.(fec)){ 
      # number of iterations
      niter <- nsamples*thin
      # print progress information at the following iteration numbers if time.start is given
      if(!verbose | is.null(time.start)){
	 printMsg <- function(...) {NULL}
      } else printMsg <- statusMsg
      niter.1percent <- c(-nburnin,floor(quantile((-nburnin):niter,prob=1:100/100)))

      ## create empty results 
      mui.samples  <- array(dim=c(nsamples, n))  
      mu.samples <- phi.samples <- numeric(nsamples)
      accept <- c("mu"=0,"phi"=0)  
      ## set initial Values
      mu <-  mu.samples[1] <- start$mu
      mui <- mui.samples[1,] <- start$mui
      phi <- phi.samples[1] <- start$phi 

      sum_mui <- sum(mui.samples[1,])
      sumlog_mui <- sum(log(mui.samples[1,]))

      # update parameters
      for (i in (-nburnin):niter){
	 # update mean
	 mu <- update_mu(mu, n*phi, sum_mui*phi)   
	 
	 # update true latent means
	 mui <- update_mui(s=fec+phi, r=prop+phi/mu)
	 sum_mui <- sum(mui)
	 sumlog_mui <- sum(log(mui))

	 # update overdispersion
	 phi <- update_phi(phi, mu, sum_mui, sumlog_mui,v.phi=v.phi)

	 # save samples
	 if(i>0 && i%%thin ==0){
	    idx <- i/thin
	    mui.samples[idx,] <- mui
	    mu.samples[idx] <- mu
	    phi.samples[idx] <- phi
	 }

	 # save acceptance
	 accept <- accept + c(attributes(mu)$accept,attributes(phi)$accept)

	 # print status information
	 if(length(w <- which(niter.1percent %in% i))) {printMsg(w-1, time.start)}
      }
      # save last value of chain
      last <- list(mu=mu, mui=mui,phi=phi)
      # return results
      res <- list(mu=mu.samples, mui=mui.samples, phi=phi.samples, accept=accept/(nburnin+1+niter),last=last)
      return(res)
   }
   )
   ####

   result <- list(mui=update_mui, mu=update_mu, phi=update_phi, runMCMC=runMCMC) 
   return(result)
}

#######################################################
### parameter updates for model PoGa2 (1-sample, zero-inflated)
#######################################################

##' Specify parameter updates for PoGa2 (1-sample) model
##'
##' Specify parameter updates for PoGa2 (1-sample) model
##'
##' @param priors  list with priors for \eqn{\mu}{mu}, \eqn{\phi}{phi} and \eqn{\psi}{psi}
##' @param fec vector with faecal egg counts
##' @param f correction factor (e.g. 50 for the McMaster technique)

##' @return a named list with update functions
##'
##' @keywords internal
setUpdates_ZIPoGa2 <- function(priors, fec, f){ 

   n <- length(fec)
   # which observed counts are zero
   whichZero_fec <- fec==0
   nNonZero_fec <- sum(fec > 0)

   #hyperparameters for mean mu
   a.mu <- priors$mu$hyperpars[1]
   b.mu <- priors$mu$hyperpars[2]
   #hyperparameters  for overdispersion parameter phi
   a.phi <- priors$phi$hyperpars[1]
   b.phi <- priors$phi$hyperpars[2]
   v.phi <- priors$phi$v
   # hyperparameters for zero-inflation parameter psi
   a.psi <- priors$psi$hyperpars[1]
   b.psi <- priors$psi$hyperpars[2]

   
   ##########################################
   ## set up update functions
   ##########################################
   # insert data/hyperparameters as default arguments
  
  
   ##################################################
   ## population mean mu
   ####################################################
   update_mu <- switch(priors$mu$proposalDist,
			"gamma" = update_h_gamma,
			"invGamma" = update_h_invGamma,
			"logNormal" = update_h_logNormal,
			"kl" = update_h_kl
			) 
   formals(update_mu)$a <- a.mu
   formals(update_mu)$b <- b.mu

   
   #############################################
   ## overdispersion parameter phi
   ##############################################
   update_phi <- switch(priors$phi$proposalDist,
			"gamma" = update_phi_gamma,
			"unif" = update_phi_unif
			)      
   formals(update_phi)$a  <- a.phi
   formals(update_phi)$b <- b.phi

   #######################################################
   ## latent means mu_i
   #####################################################
   update_mui <- update_theta_gammaMix
   formals(update_mui)$n <- n
   formals(update_mui)$whichZero_fec <- whichZero_fec
   formals(update_mui)$nNonZero_fec <- nNonZero_fec

   #######################################################
   ## zero-inflation parameter
   #####################################################
   update_psi <- update_theta_beta
   formals(update_psi)$n <- n
   formals(update_psi)$a <- a.psi
   formals(update_psi)$b <- b.psi

   ######
   runMCMC <- bquote(function(start, nsamples,nburnin, thin, v.phi,time.start=NULL, verbose=TRUE, n=.(n), nZero_fec=.(sum(fec==0)), fec=.(fec), prop=.(1/f)){
      # number of iterations
      niter <- nsamples*thin
     # print progress information at the following iteration numbers if time.start is given
      if(!verbose | is.null(time.start)){
	 printMsg <- function(...) {NULL}
      } else printMsg <- statusMsg
      niter.1percent <- c(-nburnin,floor(quantile((-nburnin):niter,prob=1:100/100)))

      ## create empty results 
      mui.samples  <- array(dim=c(nsamples, n))
      mu.samples <- phi.samples <- psi.samples <- numeric(nsamples)
      accept <- numeric(2+nZero_fec)
      names(accept) <- c("mu","phi",paste("mui",1:nZero_fec,sep=".")) 
      ## set initial Values
      mu <-  mu.samples[1] <- start$mu
      mui <- mui.samples[1,] <- start$mui
      phi <- phi.samples[1] <- start$phi 
      psi <- psi.samples[1] <- start$psi

      nNonZero_mui <- sum(mui.samples[1,]>0)
      whichNonZero_mui <- mui.samples[1,] > 0
      sum_mui <- sum(mui.samples[1,])
      sumlog_mui <- sum(log(mui.samples[1,whichNonZero_mui]))

      # update parameters
      for (i in (-nburnin):niter){
	 # update mean
	 mu <- update_mu(mu, nNonZero_mui*phi, phi*sum_mui)   

	 # update true latent means
	 mui <- update_mui(mui,fec,mu,phi,psi, p=prop)
	 nNonZero_mui <- sum(mui>0)
	 whichNonZero_mui <- mui > 0
	 sum_mui <- sum(mui)
	 sumlog_mui <- sum(log(mui[whichNonZero_mui]))
   
	 # update overdispersion
	 phi <- update_phi(phi, mu, sum_mui, sumlog_mui,nNonZero_mui, v.phi=v.phi)

	 # update zero-inflation
	 psi <- update_psi(nNonZero_mui)

	 # save samples
	 if(i>0 && i%%thin ==0){
	    idx <- i/thin
	    mui.samples[idx,] <- mui
	    mu.samples[idx] <- mu
	    phi.samples[idx] <- phi
	    psi.samples[idx] <- psi
	 }

	 # save acceptance
	 accept <- accept + c(attributes(mu)$accept,attributes(phi)$accept, attributes(mui)$accept)

	 # print status information
	 if(length(w <- which(niter.1percent %in% i))) {printMsg(w-1, time.start)}
      }
      # save last value of chain
      last <- list(mu=mu, mui=mui,phi=phi,psi=psi)
      # return results
      res <- list(mu=mu.samples,mui=mui.samples, phi=phi.samples,psi=psi.samples, accept=accept/(nburnin+1+niter),last=last)
      return(res)
   }
   )
   ####

   result <- list(mui=update_mui, mu=update_mu, phi=update_phi, psi=update_psi, runMCMC=runMCMC) 
   return(result)
}

