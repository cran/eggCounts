#############################################################
##  Faecal Egg Count Reduction testing (2-sample situation)
#############################################################

##' @include MHsteps.R
##' @include simData.R
{}

##' Set default values for the two-sample (ZI)PoGa model formulation
##'
##' @param priors.mu named list with hyperprior specifications for \eqn{\mu}{mu}, 
##'    containing elements \code{hyperpars} and \code{proposalDist} 
##' @param priors.phi named list with hyperprior specifications for the overdispersion parameter \eqn{\phi}{phi}, 
##'    containing elements \code{hyperpars}, \code{proposalDist} (unif or gamma) and 
##'    in the case of a uniform proposal also a tuning parameter \code{v} 
##' @param priors.delta named list with hyperprior specifications for the reduction in mean \eqn{\delta}{delta}, 
##'    containing elements \code{hyperpars} and \code{proposalDist}  
##' @param priors.psiB named list with hyperprior specifications for the prevalence \eqn{\psi_B}{psi.B} before treatment, 
##'    containing elements \code{hyperpars} 
##' @param priors.psiA named list with hyperprior specifications for the prevalence \eqn{\psi_A}{psi.A} after treatment, 
##'    containing elements \code{hyperpars} 
##' @param priors.deltaPsi named list with hyperprior specifications for the reduction in prevalence \eqn{\delta_\psi}{delta_psi}, 
##'    containing elements \code{hyperpars} and \code{proposalDist}  

##' @return A named list with prior specifications for \eqn{\mu}{mu}, \eqn{\phi}{phi}, \eqn{\delta}{delta} and 
##'    \eqn{\psi_B}{psi.B}, \eqn{\psi_A}{psi.A}.
##'    Default prior distributions are:
##'    \item{mu}{gamma[1, 0.001]}
##'    \item{phi}{gamma[1, 0.1]}
##'    \item{delta}{gamma[1, 1]}
##'    \item{psiB}{beta[1, 1]}
##'    \item{psiA}{beta[1, 1]}
##'    \item{deltaPsi}{beta[1, 1]}
##'   Default proposal distributions for \eqn{\mu}{mu}, \eqn{\phi}{phi} and 
##'       \eqn{\delta}{delta} (in the unpaired situation) are:
##'    \item{mu}{approximating inverse gamma, gamma or log-normal distribution}
##'    \item{phi}{unif[max(current.value - v,0), current.value + v], with \eqn{v=0.5}}
##'    \item{delta}{approximating inverse gamma, gamma or log-normal distribution}
##'
##' @keywords internal
setDefaults2 <- function(priors.mu, priors.phi, priors.delta, priors.psiB, priors.psiA, priors.deltaPsi){
   if(missing(priors.mu)) priors.mu = list(hyperpars=c(1,0.001), proposalDist="kl")
   if(is.null(priors.mu[["hyperpars", exact = TRUE]])) priors.mu$hyperpars = c(1,0.001)
   if(is.null(priors.mu[["proposalDist", exact = TRUE]])) priors.mu$proposalDist = "kl"
   
   if(missing(priors.phi)) priors.phi = list(hyperpars=c(1,0.1), proposalDist="unif", v=0.5)
   if(is.null(priors.phi[["hyperpars", exact = TRUE]])) priors.phi$hyperpars = c(1,0.1)
   if(is.null(priors.phi[["proposalDist", exact = TRUE]])) priors.phi$proposalDist = "unif"
   if(is.null(priors.phi[["v", exact = TRUE]])) priors.phi$v = 0.5

   # delta (reduction in mean) is used for PoGa model
   if(missing(priors.delta)) priors.delta <- list(priorDist="gamma", hyperpars=c(1,1), proposalDist="kl")
   if(is.null(priors.delta[["priorDist", exact = TRUE]])) priors.delta$priorDist = "gamma"
   if(is.null(priors.delta[["hyperpars", exact = TRUE]])) priors.delta$hyperpars = c(1,1)
   if(is.null(priors.delta[["proposalDist", exact = TRUE]])) priors.delta$proposalDist = "kl"

   if(missing(priors.psiB)) priors.psiB = list(hyperpars=c(1,1))
   if(is.null(priors.psiB[["hyperpars", exact = TRUE]])) priors.psiB$hyperpars = c(1,1)

   if(missing(priors.psiA)) priors.psiA = list(hyperpars=c(1,1))
   if(is.null(priors.psiA[["hyperpars", exact = TRUE]])) priors.psiA$hyperpars = c(1,1)
   
   # delta (reduction in prevalence) is used for unpaired ZIPoGa model
   if(missing(priors.deltaPsi)) priors.deltaPsi <- list(hyperpars=c(1,1), proposalDist="beta")
   if(is.null(priors.deltaPsi[["hyperpars", exact = TRUE]])) priors.deltaPsi$hyperpars = c(1,1)
   if(is.null(priors.deltaPsi[["proposalDist", exact = TRUE]])) priors.deltaPsi$proposalDist = "beta"

   return(list(mu=priors.mu,phi=priors.phi,delta=priors.delta, psiB = priors.psiB, psiA = priors.psiA, deltaPsi=priors.deltaPsi))
}

##' Set initial values for all parameters of the two-sample (ZI)PoGa model
##'
##' @param initials named list with initial values for parameters 
##'       \code{mu}  
##'       \code{phi}  
##'       \code{muiPre}  
##'       \code{muiPost}  
##'       \code{yPre}  
##'       \code{yPost}  
##'       \code{delta} 
##'       \code{psiB}  
##'       \code{psiA}  
##'       and \code{deltaPsi}  
##' @param fec.pre faecal egg count data before treatment
##' @param fec.post faecal egg count data after treatment
##' @param f.pre correction factors before treatment
##' @param f.post correction factors after treatment

##' @return A named list with initial values

##' @keywords internal
setInitials2 <- function(initials, fec.pre, fec.post, f.pre, f.post){
   # set initials
   yPre <- fec.pre*f.pre
   yPost <- fec.post*f.post
   mu <-  mean(yPre)
   muiPre <- yPre + rbinom(length(fec.pre), size = 1, prob = 0.5)
   muiPost <- yPost + rbinom(length(fec.post), size = 1, prob = 0.5)
   phi <- .8 
   delta <- .5
   psiB <- max(0.1, mean(fec.pre==1))
   psiA <- max(0.1, mean(fec.post==1))
   deltaPsi <- .5

   if(is.null(initials)) initials <- list(mu=mu, phi=phi, yPre=yPre, yPost=yPost, muiPre=muiPre, muiPost=muiPost, delta=delta, psiB=psiB, psiA=psiA)

   if(is.null(initials[["yPre", exact = TRUE]])) initials$yPre  <- yPre
   if(is.null(initials[["yPost", exact = TRUE]])) initials$yPost  <- yPost
   if(is.null(initials[["mu", exact = TRUE]])) initials$mu  <- mu
   if(is.null(initials[["muiPre", exact = TRUE]])) initials$muiPre  <- muiPre
   if(is.null(initials[["muiPost", exact = TRUE]])) initials$muiPost  <- muiPost
   if(is.null(initials[["phi", exact = TRUE]])) initials$phi  <- phi
   if(is.null(initials[["delta", exact = TRUE]])) initials$delta  <- delta
   if(is.null(initials[["psiB", exact = TRUE]])) initials$psiB  <- psiB
   if(is.null(initials[["psiA", exact = TRUE]])) initials$psiA  <- psiA
   if(is.null(initials[["deltaPsi", exact = TRUE]])) initials$deltaPsi  <- deltaPsi

   # ToDo: check user-defined initials
   
   return(initials)
}



###################################
#  PoGa (paired and unpaired)

##' Modelling of Faecal Egg Count data (two-sample case)
##'
##' Modelling of Faecal Egg Count data in a two-sample case using
##' a paired or unpaired (ZI) Poisson gamma model formulation
##'
##' @param fec.pre vector with faecal egg counts before treatment
##' @param fec.post vector with faecal egg counts after treatment
##' @param rawCounts logical indicating whether \code{fec.pre} and \code{fec.post} corresponds to raw counts
##'  (as counted on the McMaster slide), or to calculated EpGs (raw counts times correction factor). 
##'  Defaults to \code{FALSE}.
##' @param f.pre correction factor(s) before treatment
##' @param f.post correction factor(s) after treatment 
##' @param model string with model formulation ("paired" or "unpaired")
##' @param priors.mu list with hyper-prior/proposal information for \eqn{\mu}{mu}
##' @param priors.phi list with hyper-prior/proposal information for \eqn{\phi}{phi}
##' @param priors.delta list with hyper-prior/proposal information for \eqn{\delta}{delta}
##' @param priors.psiA list with hyper-prior/proposal information for \eqn{\psi_A}{psiA}
##' @param priors.psiB list with hyper-prior/proposal information for \eqn{\psi_B}{psiB}
##' @param priors.deltaPsi list with hyper-prior/proposal information for \eqn{\delta_\psi}{deltaPsi}
##' @param maxiter.pilot maximal number of tries to determine a good tuning value
##'   for the proposal distribution for \eqn{\phi}{phi}
##' @param nburnin number of burn-in iterations
##' @param nsamples number of desired samples
##' @param thin thinning parameter
##' @param initials named list with starting values for the parameters 
##'         \code{mu}, \code{phi}, \code{muiPre}, \code{muiPost}, \code{yPre}, \code{yPost},
##'         \code{delta}, \code{psiB} and \code{psiA}
##' @param verbose print progress information
##' @param \dots extra arguments (not used atm)

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
# @TODO testing of a zero-inflated model

##' @export
fecr_mcmc <- function(
   fec.pre,    # FEC data pre-treatment
   fec.post,   # FEC data post-treatment
   rawCounts = FALSE,
   f.pre = 50, # correction factor 
   f.post = f.pre,
   model = "paired",  # assumed model: "unpaired","ZIPoGa_u","ZIPoGa_upd"
   priors.mu = list(hyperpars= c(1,0.001),
		    proposalDist = "kl"),
   priors.phi = list(hyperpars = c(1,0.1),
		     proposalDist = "unif",
		     v = 0.5   # tuning param for uniform proposal
		    ),
   priors.delta = list(priorDist = c("gamma","beta")[1],
		       hyperpars = c(1,1),
	               proposalDist = NULL),
   priors.psiB = list(hyperpars= c(1,1)),
   priors.psiA = list(hyperpars= c(1,1)),
   priors.deltaPsi = list(hyperpars= c(1,1), proposalDist = "beta"),
   maxiter.pilot = 15,  # number of tuning tries
   nburnin = 5e3,
   nsamples = 1e4,
   thin =1,
   initials=NULL,
   verbose=TRUE,
   ...
   ){

	## TODO: check data format

   # number of faecal samples
   n.pre <- length(fec.pre)
   n.post <- length(fec.post)

   # check correction factors
   if(length(f.pre)>1 && length(f.pre)!=n.pre) stop("Lengths of the vectors f.pre and fec.pre do not match\n")
   if(length(f.post)>1 && length(f.post)!=n.post) stop("Lengths of the vectors f.post and fec.post do not match\n")

   # raw counts or EpGs?
   dilution.pre <- f.pre; dilution.post <- f.post
   if(rawCounts){
      dilution.pre <- dilution.post <- 1
   } 
   # divide data by correction factor
   if(any(fec.pre %% dilution.pre !=0)) stop("Correction factor f.pre does not match the given pre-treatment fec\n")
   fec.pre <- fec.pre/dilution.pre
   if(any(fec.post %% dilution.post !=0)) stop("Correction factor f.post does not match the given post-treatment fec\n")
   fec.post <- fec.post/dilution.post


   # check model and set default values
   model <- match.arg(model, c("paired","unpaired","ZIPoGa_u","ZIPoGa_upd"))
   priors <- setDefaults2(priors.mu=priors.mu, priors.phi=priors.phi, priors.delta=priors.delta, priors.psiB=priors.psiB, priors.psiA=priors.psiA, priors.deltaPsi=priors.deltaPsi)

   # set update functions
   updates <- switch(model,
                    "paired"= setUpdates_PoGa_p(priors, fec.pre, fec.post, f.pre, f.post),
                    "unpaired"= setUpdates_PoGa_u(priors, fec.pre, fec.post, f.pre, f.post),
                    "ZIPoGa_u"= setUpdates_ZIPoGa_u(priors, fec.pre, fec.post, f.pre, f.post),
                    "ZIPoGa_upd"= setUpdates_ZIPoGa_u_pd(priors, fec.pre, fec.post, f.pre, f.post)
              )
   if(model %in% c("paired") && n.pre != n.post){
     stop("post sample size different to pre sample size\n")
   }

   update_yPre <- updates$yPre
   update_yPost <- updates$yPost
   update_muiPre <- updates$muiPre
   update_muiPost <- updates$muiPost
   update_mu <- updates$mu
   update_delta <- updates$delta
   update_phi <- updates$phi
   update_psiB <- updates$psiB
   update_psiA <- updates$psiA
   update_deltaPsi <- updates$deltaPsi
   runMCMC <- eval(updates$runMCMC)
   
   # set initials 
   start <- setInitials2(initials, fec.pre, fec.post, f.pre, f.post)
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
   # update parameters
   time.start <- Sys.time()
   samples <- runMCMC(start,nsamples=nsamples,nburnin=nburnin,thin=thin, v.phi=v.phi, time.start=time.start, verbose=verbose)

   #return results
   result <- list(samples=samples, priors=priors, v.phi=v.phi,initials=start, model=model, nburnin=nburnin, nsamples=nsamples,thin=thin)

   class(result) <- "fecrm"
   return(result)
}


######
### parameter updates for model PoGa_p: paired sample situation
#######

##' Specify parameter updates for PoGa_p (paired 2-sample) model
##'
##' @param priors  list with prior specifications
##' @param fec.pre vector with faecal egg counts (before treatment)
##' @param fec.post vector with faecal egg counts (after treatment)
##' @param f.pre correction factors before treatment 
##' @param f.post correction factors after treatment 

##' @return a named list with update functions
##'
##' @keywords internal
setUpdates_PoGa_p <- function(priors, fec.pre, fec.post, f.pre, f.post){ 

   n <- length(fec.pre)
   propC.pre <- 1- 1/f.pre
   propC.post <- 1- 1/f.post

   #hyperparameters for pre-treatment mean mu
   a.mu <- priors$mu$hyperpars[1]
   b.mu <- priors$mu$hyperpars[2]
   #hyperparameters  for overdispersion parameter phi
   a.phi <- priors$phi$hyperpars[1]
   b.phi <- priors$phi$hyperpars[2]
   v.phi <- priors$phi$v
   #hyperparameters  for change in mean delta
   a.delta <- priors$delta$hyperpars[1]
   b.delta <- priors$delta$hyperpars[2]

#check that (a,b)-values are positive
   

   ## set up update functions
   ## insert data/hyperparameters as default arguments

   ### latent true egg counts y_i pre-treamtent
   update_yPre <- update_Y
   formals(update_yPre)$n <- n
   formals(update_yPre)$propC <- propC.pre
   formals(update_yPre)$fec <- fec.pre

   ### latent true egg counts y_i post-treamtent
   update_yPost <- update_Y
   formals(update_yPost)$n <- n
   formals(update_yPost)$propC <- propC.post
   formals(update_yPost)$fec <- fec.post

   ### pre-treatment mean mu
   update_mu <- switch(priors$mu$proposalDist,
			"gamma" = update_h_gamma,
			"invGamma" = update_h_invGamma,
			"logNormal" = update_h_logNormal,
			"kl" = update_h_kl
			)    
   formals(update_mu)$a <- a.mu
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

   ### mean epg rates before treatment
   update_muiPre <- update_theta_gamma
   formals(update_muiPre)$n <- n

   ### reduction in mean
   #update_delta <- update_theta_gamma1
   update_delta <- switch(priors$delta$priorDist,
			"beta" = update_delta_beta,
			"gamma" = update_theta_gamma1         # in case of gamma prior, FC is gamma
			) 
   formals(update_delta)$a <- a.delta
   formals(update_delta)$b <- b.delta

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
      yPre.samples <-yPost.samples <- muiPre.samples <- array(dim=c(nsamples, n))
      mu.samples <- phi.samples <- delta.samples <- numeric(nsamples)
      accept <- numeric(2+ (priors$delta$priorDist=="beta"))
      names(accept) <- c("mu","phi","delta")[1:length(accept)] 
 
      ## set initials values
      yPre <- yPre.samples[1,] <- start$yPre
      yPost <- yPost.samples[1,] <- start$yPost
      mu <-  mu.samples[1] <- start$mu
      delta <- delta.samples[1] <- start$delta
      muiPre <- muiPre.samples[1,] <- start$muiPre
      phi <- phi.samples[1] <- start$phi

      ## compute frequently used quantities for updates
      sum_muiPre <- sum(muiPre.samples[1,])
      sumlog_muiPre <- sum(log(muiPre.samples[1,]))

      # update parameters
      for (i in (-nburnin):niter){
	 # update mean
	 mu <- update_mu(mu, n*phi, phi*sum_muiPre)  
	 
	 # update latent epg counts
	 yPre <- update_yPre(muiPre)
	 yPost <- update_yPost(muiPre*delta)

	 # update true latent means
	 muiPre <- update_muiPre(s=yPre+yPost+phi, r=1+delta+phi/mu) 
	 sum_muiPre <- sum(muiPre)
	 sumlog_muiPre <- sum(log(muiPre))

	 # update overdispersion
	 phi <- update_phi(phi, mu, sum_muiPre, sumlog_muiPre,v.phi=v.phi)

	 # update change in mean
	 #delta <- update_delta(sum(yPost), sum_muiPre)
	 delta <- update_delta(delta=delta, sum(yPost), sum_muiPre)

	 # save samples
	 if(i>0 && i%%thin ==0){
	    idx <- i/thin
	    yPre.samples[idx,] <- yPre
	    yPost.samples[idx,] <- yPost
	    muiPre.samples[idx,] <- muiPre
	    delta.samples[idx] <- delta
	    mu.samples[idx] <- mu
	    phi.samples[idx] <- phi
	 }

	 # save acceptance
	 accept <- accept + c(attributes(mu)$accept,attributes(phi)$accept, attributes(delta)$accept)

	 # print status information
	 if(length(w <- which(niter.1percent %in% i))) {printMsg(w-1, time.start)}
      }
      # save last value of chain
      last <- list(yPre=yPre,yPost=yPost, mu=mu, muiPre=muiPre,phi=phi,delta=delta)
      # return results
      res <- list(mu=mu.samples, yPre=yPre.samples, yPost=yPost.samples, muiPre=muiPre.samples, phi=phi.samples, delta=delta.samples, accept=accept/(nburnin+1+niter),last=last)
      return(res)
   }
   )
   #####

   result <- list(yPre=update_yPre, yPost=update_yPost,muiPre=update_muiPre,mu=update_mu, delta=update_delta, phi=update_phi, runMCMC=runMCMC) 
   return(result)
}

######
### parameter updates for model PoGa_u: unpaired 2-sample situation
#######

##' Specify parameter updates for PoGa_u (unpaired 2-sample) model
##'
##' @param priors  list with prior specifications
##' @param fec.pre vector with faecal egg counts (before treatment)
##' @param fec.post vector with faecal egg counts (after treatment)
##' @param f.pre correction factors before treatment 
##' @param f.post correction factors after treatment 

##' @return a named list with update functions
##'
##' @keywords internal
setUpdates_PoGa_u <- function(priors, fec.pre, fec.post, f.pre, f.post){ 

   n.pre <- length(fec.pre)
   n.post <- length(fec.post)
   propC.pre <- 1- 1/f.pre
   propC.post <- 1- 1/f.post

   #hyperparameters for pre-treatment mean mu
   a.mu <- priors$mu$hyperpars[1]
   b.mu <- priors$mu$hyperpars[2]
   #hyperparameters  for overdispersion parameter phi
   a.phi <- priors$phi$hyperpars[1]
   b.phi <- priors$phi$hyperpars[2]
   v.phi <- priors$phi$v
   #hyperparameters  for change in mean delta
   a.delta <- priors$delta$hyperpars[1]
   b.delta <- priors$delta$hyperpars[2]

#check that (a,b)-values are positive
   

   ## set up update functions
   ## insert data/hyperparameters as default arguments

   ### latent true egg counts y_i pre-treamtent
   update_yPre <- update_Y
   formals(update_yPre)$n <- n.pre
   formals(update_yPre)$propC <- propC.pre
   formals(update_yPre)$fec <- fec.pre

   ### latent true egg counts y_i post-treamtent
   update_yPost <- update_Y
   formals(update_yPost)$n <- n.post
   formals(update_yPost)$propC <- propC.post
   formals(update_yPost)$fec <- fec.post


   ### pre-treatment mean mu
   update_muiPre <- update_theta_gamma
   formals(update_muiPre)$n <- n.pre

   ### post-treatment mean mu
   update_muiPost <- update_theta_gamma
   formals(update_muiPost)$n <- n.post

   ## change in mean
   update_delta <- switch(priors$delta$proposalDist,
			"gamma" = update_h_gamma,
			"invGamma" = update_h_invGamma,
			"logNormal" = update_h_logNormal,
			"kl" = update_h_kl
			)    
   formals(update_delta)$a <- a.delta
   formals(update_delta)$b <- b.delta

   
   ### pre-treatment mean mu
   update_mu <- switch(priors$mu$proposalDist,
			"gamma" = update_h_gamma,
			"invGamma" = update_h_invGamma,
			"logNormal" = update_h_logNormal,
			"kl" = update_h_kl
			)    
   formals(update_mu)$a <- a.mu
   formals(update_mu)$b <- b.mu


   ### overdispersion parameter phi
   update_phi <- switch(priors$phi$proposalDist,
			"gamma" = update_phi_gamma,
			"unif" = update_phi_unif,
			"unif01" = update_phiTransf_unif
			)      
   formals(update_phi)$a  <- a.phi
   formals(update_phi)$b <- b.phi
   formals(update_phi)$n <- n.pre+n.post
     

   ######
   runMCMC <- bquote(function(start, nsamples,nburnin,thin, v.phi,time.start=NULL, verbose=TRUE, n.pre=.(n.pre), n.post=.(n.post)){ 
      # total sample size FEC
      n.all <- n.pre+n.post
      
      # number of iterations
      niter <- nsamples*thin
      # print progress information at the following iteration numbers if time.start is given
      if(!verbose | is.null(time.start)){
	 printMsg <- function(...) {NULL}
      } else printMsg <- statusMsg
      niter.1percent <- c(-nburnin,floor(quantile((-nburnin):niter,prob=1:100/100)))

      ## create empty results 
      yPre.samples <- muiPre.samples <- array(dim=c(nsamples, n.pre))
      yPost.samples <- muiPost.samples <- array(dim=c(nsamples, n.post))
      mu.samples <- phi.samples <- delta.samples <- numeric(nsamples)
      accept <- numeric(3)
      names(accept) <- c("mu","phi","delta") 
 
      ## set initials values
      yPre <- yPre.samples[1,] <- start$yPre
      yPost <- yPost.samples[1,] <- start$yPost
      mu <-  mu.samples[1] <- start$mu
      delta <- delta.samples[1] <- start$delta
      muiPre <- muiPre.samples[1,] <- start$muiPre
      muiPost <- muiPost.samples[1,] <- start$muiPost
      phi <- phi.samples[1] <- start$phi

      ## compute frequently used quantities for updates
      sum_muiPre <- sum(muiPre.samples[1,])
      sumlog_muiPre <- sum(log(muiPre.samples[1,]))
      sum_muiPost <- sum(muiPost.samples[1,])
      sumlog_muiPost <- sum(log(muiPost.samples[1,]))

      # update parameters
      for (i in (-nburnin):niter){
	 # update mean
	 mu <- update_mu(mu, n.all*phi, phi*(sum_muiPre+sum_muiPost/delta))  
	 
	 # update latent epg counts
	 yPre <- update_yPre(muiPre)
	 yPost <- update_yPost(muiPost)

	 # update true latent means
	 muiPre <- update_muiPre(s=yPre+phi, r=1+phi/mu)
	 sum_muiPre <- sum(muiPre)
	 sumlog_muiPre <- sum(log(muiPre))

	 muiPost <- update_muiPost(s=yPost+phi, r=1+phi/mu/delta)
	 sum_muiPost <- sum(muiPost)
	 sumlog_muiPost <- sum(log(muiPost))

	 # update overdispersion
	 phi <- update_phi(phi, mu, sum_mui=sum_muiPre+sum_muiPost/delta, sumlog_mui=sumlog_muiPre+sumlog_muiPost, v.phi=v.phi)

	 # update change in mean
	 delta <- update_delta(delta, n.post*phi, sum_muiPost*phi/mu) 

	 # save samples
	 if(i>0 && i%%thin ==0){
	    idx <- i/thin
	    yPre.samples[idx,] <- yPre
	    yPost.samples[idx,] <- yPost
	    muiPre.samples[idx,] <- muiPre
	    muiPost.samples[idx,] <- muiPost
	    delta.samples[idx] <- delta
	    mu.samples[idx] <- mu
	    phi.samples[idx] <- phi
	 }

	 # save acceptance
	 accept <- accept + c(attributes(mu)$accept,attributes(phi)$accept,attributes(delta)$accept)

	 # print status information
	 if(length(w <- which(niter.1percent %in% i))) {printMsg(w-1, time.start)}
      }
      # save last value of chain
      last <- list(yPre=yPre,yPost=yPost, mu=mu, muiPre=muiPre, muiPost=muiPost,phi=phi,delta=delta)
      # return results
      res <- list(mu=mu.samples, yPre=yPre.samples, yPost=yPost.samples, muiPre=muiPre.samples, muiPost=muiPost.samples, phi=phi.samples, delta=delta.samples, accept=accept/(nburnin+1+niter),last=last)
      return(res)
   }
   )
   #####


   result <- list(yPre=update_yPre, yPost=update_yPost,muiPre=update_muiPre,muiPost=update_muiPost, mu=update_mu, delta=update_delta, phi=update_phi, runMCMC=runMCMC) 
   return(result)
}

###########################################################
### parameter updates for model ZIPoGa_u: unpaired sample situation, reduced prevalence
###########################################################

##' Specify parameter updates for ZIPoGa_u (2-sample) model
##'
##' @param priors  list with prior specifications
##' @param fec.pre vector with faecal egg counts (before treatment)
##' @param fec.post vector with faecal egg counts (after treatment)
##' @param f.pre correction factors before treatment 
##' @param f.post correction factors after treatment 

##' @return a named list with update functions
##'
##' @keywords internal
setUpdates_ZIPoGa_u <- function(priors, fec.pre, fec.post, f.pre, f.post){ 

   n.pre <- length(fec.pre)
   n.post <- length(fec.post)
   # which observed counts are zero
   whichZero_fec.pre <- fec.pre==0
   nNonZero_fec.pre <- sum(fec.pre > 0)
   nZero_fec.pre <- sum(fec.pre==0)
   whichZero_fec.post <- fec.post==0
   nNonZero_fec.post <- sum(fec.post > 0)
   nZero_fec.post <- sum(fec.post==0)

   propC.pre <- 1- 1/f.pre
   propC.post <- 1- 1/f.post

   #hyperparameters for pre-treatment mean mu
   a.mu <- priors$mu$hyperpars[1]
   b.mu <- priors$mu$hyperpars[2]
   #hyperparameters  for overdispersion parameter phi
   a.phi <- priors$phi$hyperpars[1]
   b.phi <- priors$phi$hyperpars[2]
   v.phi <- priors$phi$v
   # hyperparameters for zero-inflation parameter psi
   a.psiB <- priors$psiB$hyperpars[1]
   b.psiB <- priors$psiB$hyperpars[2]
   a.psiA <- priors$psiA$hyperpars[1]
   b.psiA <- priors$psiA$hyperpars[2]

#check that (a,b)-values are positive
   

   ## set up update functions
   ## insert data/hyperparameters as default arguments

   ### latent true egg counts y_i pre-treamtent
   update_yPre <- update_Y
   formals(update_yPre)$n <- n.pre
   formals(update_yPre)$propC <- propC.pre
   formals(update_yPre)$fec <- fec.pre

   ### latent true egg counts y_i post-treamtent
   update_yPost <- update_Y
   formals(update_yPost)$n <- n.post
   formals(update_yPost)$propC <- propC.post
   formals(update_yPost)$fec <- fec.post


   ### pre-treatment mean mu
   update_mu <- switch(priors$mu$proposalDist,
			"gamma" = update_h_gamma,
			"invGamma" = update_h_invGamma,
			"logNormal" = update_h_logNormal,
			"kl" = update_h_kl
			)    
   formals(update_mu)$a <- a.mu
   formals(update_mu)$b <- b.mu

   ### overdispersion parameter phi
   update_phi <- switch(priors$phi$proposalDist,
			"gamma" = update_phi_gamma,
			"unif" = update_phi_unif,
			"unif01" = update_phiTransf_unif
			)      
   formals(update_phi)$a  <- a.phi
   formals(update_phi)$b <- b.phi


   ### pre-treatment mean epg rates
   update_muiPre <- update_theta_gammaMix
   formals(update_muiPre)$n <- n.pre
   formals(update_muiPre)$whichZero_fec <- whichZero_fec.pre
   formals(update_muiPre)$nNonZero_fec <- nNonZero_fec.pre

   ### post-treatment mean epg rates
   update_muiPost <- update_theta_gammaMix
   formals(update_muiPost)$n <- n.post
   formals(update_muiPost)$whichZero_fec <- whichZero_fec.post
   formals(update_muiPost)$nNonZero_fec <- nNonZero_fec.post
 
   ### zero-inflation parameter
   update_psiB <- update_theta_beta
   formals(update_psiB)$a <- a.psiB
   formals(update_psiB)$b <- b.psiB
   formals(update_psiB)$n <- n.pre

   update_psiA <- update_theta_beta
   formals(update_psiA)$a <- a.psiA
   formals(update_psiA)$b <- b.psiA
   formals(update_psiA)$n <- n.post
    

   ######
   runMCMC <- bquote(function(start, nsamples,nburnin,thin, v.phi,time.start=NULL, verbose=TRUE, n.pre=.(n.pre), n.post=.(n.post), nZero_fec.pre=.(nZero_fec.pre), nZero_fec.post=.(nZero_fec.post)){ 
      # total sample size FEC
      n.all <- n.pre+n.post
      # number of iterations
      niter <- nsamples*thin
      # print progress information at the following iteration numbers if time.start is given
      if(!verbose | is.null(time.start)){
	 printMsg <- function(...) {NULL}
      } else printMsg <- statusMsg
      niter.1percent <- c(-nburnin,floor(quantile((-nburnin):niter,prob=1:100/100)))

      ## create empty results 
      yPre.samples <- muiPre.samples <- array(dim=c(nsamples, n.pre))
      yPost.samples <- muiPost.samples <- array(dim=c(nsamples, n.post))
      mu.samples <- phi.samples <- psiB.samples <-psiA.samples <- numeric(nsamples)
      accept <- numeric(2+nZero_fec.pre+nZero_fec.post)
      namsAccept <- c("mu","phi")
      if(nZero_fec.pre>0) namsAccept <- c(namsAccept, paste("muiB", 1:nZero_fec.pre,sep="."))
      if(nZero_fec.post>0) namsAccept <- c(namsAccept, paste("muiA", 1:nZero_fec.post,sep="."))
      names(accept) <- namsAccept #c("mu","phi",paste("muiB", 1:nZero_fec.pre,sep="."), paste("muiA", 1:nZero_fec.post,sep=".")) 
 
      ## set initials values
      yPre <- yPre.samples[1,] <- start$yPre
      yPost <- yPost.samples[1,] <- start$yPost
      mu <-  mu.samples[1] <- start$mu
      muiPre <- muiPre.samples[1,] <- start$muiPre
      muiPost <- muiPost.samples[1,] <- start$muiPost
      phi <- phi.samples[1] <- start$phi
      psiB <- psiB.samples[1] <- start$psiB
      psiA <- psiA.samples[1] <- start$psiA

      ## compute frequently used quantities for updates
      nNonZero_muiPre <- sum(muiPre.samples[1,]>0)
      nNonZero_muiPost <- sum(muiPost.samples[1,]>0)
      nNonZero_muiAll <- nNonZero_muiPre + nNonZero_muiPost

      whichNonZero_muiPre <- muiPre.samples[1,] > 0
      whichNonZero_muiPost <- muiPost.samples[1,] > 0

      sum_muiPre <- sum(muiPre.samples[1,])
      sum_muiPost <- sum(muiPost.samples[1,])
      sum_muiAll <- sum_muiPre+sum_muiPost

      sumlog_muiPre <- sum(log(muiPre.samples[1,whichNonZero_muiPre]))
      sumlog_muiPost <- sum(log(muiPost.samples[1,whichNonZero_muiPost]))
      sumlog_muiAll <- sumlog_muiPre+sum_muiPost

      # update parameters
      for (i in (-nburnin):niter){
	 # update mean
	 mu <- update_mu(mu, nNonZero_muiAll*phi, phi*sum_muiAll)  
	 
	 # update latent epg counts
	 yPre <- update_yPre(muiPre)
	 yPost <- update_yPost(muiPost)

	 # update true latent means
	 muiPre <- update_muiPre(muiPre,yPre,mu,phi,psiB)
	 nNonZero_muiPre <- sum(muiPre>0)
	 whichNonZero_muiPre <- muiPre > 0
	 sum_muiPre <- sum(muiPre)
	 sumlog_muiPre <- sum(log(muiPre[whichNonZero_muiPre]))

	 muiPost <- update_muiPost(muiPost,yPost,mu,phi,psiA)
	 nNonZero_muiPost <- sum(muiPost>0)
	 whichNonZero_muiPost <- muiPost > 0
	 sum_muiPost <- sum(muiPost)
	 sumlog_muiPost <- sum(log(muiPost[whichNonZero_muiPost]))

	 nNonZero_muiAll <- nNonZero_muiPre +nNonZero_muiPost
	 sum_muiAll <- sum_muiPre+sum_muiPost
	 sumlog_muiAll <- sumlog_muiPre+sumlog_muiPost

	 # update overdispersion
	 phi <- update_phi(phi, mu, sum_muiAll, sumlog_muiAll, nNonZero_muiAll,v.phi=v.phi)

	 # update zero-inflation
	 psiB <- update_psiB(nNonZero_muiPre)
	 psiA <- update_psiA(nNonZero_muiPost)

	 # save samples
	 if(i>0 && i%%thin ==0){
	    idx <- i/thin
	    yPre.samples[idx,] <- yPre
	    yPost.samples[idx,] <- yPost
	    muiPre.samples[idx,] <- muiPre
	    muiPost.samples[idx,] <- muiPost
	    mu.samples[idx] <- mu
	    phi.samples[idx] <- phi
	    psiB.samples[idx] <- psiB
	    psiA.samples[idx] <- psiA
	 }

	 # save acceptance
	 accept <- accept + c(attributes(mu)$accept,attributes(phi)$accept, attributes(muiPre)$accept, attributes(muiPost)$accept) 

	 # print status information
	 if(length(w <- which(niter.1percent %in% i))) {printMsg(w-1, time.start)}
      }
      # save last value of chain
      last <- list(yPre=yPre,yPost=yPost, mu=mu, muiPre=muiPre, muiPost=muiPost,phi=phi,psiB=psiB, psiA=psiA)
      # return results
      res <- list(mu=mu.samples, yPre=yPre.samples, yPost=yPost.samples, muiPre=muiPre.samples, muiPost=muiPost.samples, phi=phi.samples, psiPre=psiB.samples, psiPost=psiA.samples, accept=accept/(nburnin+1+niter),last=last)
      return(res)
   }
   )
   #####


   result <- list(yPre=update_yPre, yPost=update_yPost,muiPre=update_muiPre,muiPost=update_muiPost, mu=update_mu, psiB=update_psiB, psiA=update_psiA, phi=update_phi, runMCMC=runMCMC) 
   return(result)
}



###########################################################
### parameter updates for model ZIPoGa_u: unpaired sample situation, reduced prevalence
###########################################################

##' Specify parameter updates for ZIPoGa_u_pd (2-sample) model
##'
##' @param priors  list with prior specifications
##' @param fec.pre vector with faecal egg counts (before treatment)
##' @param fec.post vector with faecal egg counts (after treatment)
##' @param f.pre correction factors before treatment 
##' @param f.post correction factors after treatment 

##' @return a named list with update functions
##'
##' @keywords internal
setUpdates_ZIPoGa_u_pd <- function(priors, fec.pre, fec.post, f.pre, f.post){ 

   n.pre <- length(fec.pre)
   n.post <- length(fec.post)
   # which observed counts are zero
   whichZero_fec.pre <- fec.pre==0
   nNonZero_fec.pre <- sum(fec.pre > 0)
   nZero_fec.pre <- sum(fec.pre==0)
   whichZero_fec.post <- fec.post==0
   nNonZero_fec.post <- sum(fec.post > 0)
   nZero_fec.post <- sum(fec.post==0)

   propC.pre <- 1- 1/f.pre
   propC.post <- 1- 1/f.post

   #hyperparameters for pre-treatment mean mu
   a.mu <- priors$mu$hyperpars[1]
   b.mu <- priors$mu$hyperpars[2]
   #hyperparameters  for overdispersion parameter phi
   a.phi <- priors$phi$hyperpars[1]
   b.phi <- priors$phi$hyperpars[2]
   v.phi <- priors$phi$v
   # hyperparameters for zero-inflation parameter psi
   a.psiB <- priors$psiB$hyperpars[1]
   b.psiB <- priors$psiB$hyperpars[2]
   # hyperparameters for reduction in prevalence deltaPsi
   a.deltaPsi <- priors$deltaPsi$hyperpars[1]
   b.deltaPsi <- priors$deltaPsi$hyperpars[2]

#check that (a,b)-values are positive
   

   ## set up update functions
   ## insert data/hyperparameters as default arguments

   ### latent true egg counts y_i pre-treamtent
   update_yPre <- update_Y
   formals(update_yPre)$n <- n.pre
   formals(update_yPre)$propC <- propC.pre
   formals(update_yPre)$fec <- fec.pre

   ### latent true egg counts y_i post-treamtent
   update_yPost <- update_Y
   formals(update_yPost)$n <- n.post
   formals(update_yPost)$propC <- propC.post
   formals(update_yPost)$fec <- fec.post


   ### pre-treatment mean mu
   update_mu <- switch(priors$mu$proposalDist,
			"gamma" = update_h_gamma,
			"invGamma" = update_h_invGamma,
			"logNormal" = update_h_logNormal,
			"kl" = update_h_kl
			)    
   formals(update_mu)$a <- a.mu
   formals(update_mu)$b <- b.mu

   ### overdispersion parameter phi
   update_phi <- switch(priors$phi$proposalDist,
			"gamma" = update_phi_gamma,
			"unif" = update_phi_unif,
			"unif01" = update_phiTransf_unif
			)      
   formals(update_phi)$a  <- a.phi
   formals(update_phi)$b <- b.phi


   ### pre-treatment mean epg rates
   update_muiPre <- update_theta_gammaMix
   formals(update_muiPre)$n <- n.pre
   formals(update_muiPre)$whichZero_fec <- whichZero_fec.pre
   formals(update_muiPre)$nNonZero_fec <- nNonZero_fec.pre

   ### post-treatment mean epg rates
   update_muiPost <- update_theta_gammaMix
   formals(update_muiPost)$n <- n.post
   formals(update_muiPost)$whichZero_fec <- whichZero_fec.post
   formals(update_muiPost)$nNonZero_fec <- nNonZero_fec.post
 
   ### zero-inflation parameter
   update_psiB <- update_psi_beta
   formals(update_psiB)$a <- a.psiB
   formals(update_psiB)$b <- b.psiB
   formals(update_psiB)$n <- n.post

   update_deltaPsi <- update_psi_beta
   formals(update_deltaPsi)$a <- a.deltaPsi
   formals(update_deltaPsi)$b <- b.deltaPsi
   formals(update_deltaPsi)$n <- n.post
    

   ######
   runMCMC <- bquote(function(start, nsamples,nburnin,thin, v.phi,time.start=NULL, verbose=TRUE, n.pre=.(n.pre), n.post=.(n.post), nZero_fec.pre=.(nZero_fec.pre), nZero_fec.post=.(nZero_fec.post)){ 
      # total sample size FEc
      n.all <- n.pre+n.post

      # number of iterations
      niter <- nsamples*thin
      # print progress information at the following iteration numbers if time.start is given
      if(!verbose | is.null(time.start)){
	 printMsg <- function(...) {NULL}
      } else printMsg <- statusMsg
      niter.1percent <- c(-nburnin,floor(quantile((-nburnin):niter,prob=1:100/100)))

      ## create empty results 
      yPre.samples <- muiPre.samples <- array(dim=c(nsamples, n.pre))
      yPost.samples <- muiPost.samples <- array(dim=c(nsamples, n.post))
      mu.samples <- phi.samples <- psiB.samples <-deltaPsi.samples <- numeric(nsamples)
      accept <- numeric(4+nZero_fec.pre+nZero_fec.post)
      namsAccept <- c("mu","phi","psiB","deltaPsi")
      if(nZero_fec.pre>0) namsAccept <- c(namsAccept, paste("muiB", 1:nZero_fec.pre,sep="."))
      if(nZero_fec.post>0) namsAccept <- c(namsAccept, paste("muiA", 1:nZero_fec.post,sep="."))
      names(accept) <- namsAccept #c("mu","phi","psiB","deltaPsi", paste("muiB", 1:nZero_fec.pre,sep="."), paste("muiA", 1:nZero_fec.post,sep=".")) 
 
      ## set initials values
      yPre <- yPre.samples[1,] <- start$yPre
      yPost <- yPost.samples[1,] <- start$yPost
      mu <-  mu.samples[1] <- start$mu
      muiPre <- muiPre.samples[1,] <- start$muiPre
      muiPost <- muiPost.samples[1,] <- start$muiPost
      phi <- phi.samples[1] <- start$phi
      psiB <- psiB.samples[1] <- start$psiB
      deltaPsi <- deltaPsi.samples[1] <- start$deltaPsi

      ## compute frequently used quantities for updates
      nNonZero_muiPre <- sum(muiPre.samples[1,]>0)
      nNonZero_muiPost <- sum(muiPost.samples[1,]>0)
      nNonZero_muiAll <- nNonZero_muiPre + nNonZero_muiPost

      whichNonZero_muiPre <- muiPre.samples[1,] > 0
      whichNonZero_muiPost <- muiPost.samples[1,] > 0

      sum_muiPre <- sum(muiPre.samples[1,])
      sum_muiPost <- sum(muiPost.samples[1,])
      sum_muiAll <- sum_muiPre+sum_muiPost

      sumlog_muiPre <- sum(log(muiPre.samples[1,whichNonZero_muiPre]))
      sumlog_muiPost <- sum(log(muiPost.samples[1,whichNonZero_muiPost]))
      sumlog_muiAll <- sumlog_muiPre+sum_muiPost

      # update parameters
      for (i in (-nburnin):niter){
	 # update mean
	 mu <- update_mu(mu, nNonZero_muiAll*phi, phi*sum_muiAll)  
	 
	 # update latent epg counts
	 yPre <- update_yPre(muiPre)
	 yPost <- update_yPost(muiPost)

	 # update true latent means
	 muiPre <- update_muiPre(muiPre,yPre,mu,phi,psiB)
	 nNonZero_muiPre <- sum(muiPre>0)
	 whichNonZero_muiPre <- muiPre > 0
	 sum_muiPre <- sum(muiPre)
	 sumlog_muiPre <- sum(log(muiPre[whichNonZero_muiPre]))

	 muiPost <- update_muiPost(muiPost,yPost,mu,phi,deltaPsi*psiB)
	 nNonZero_muiPost <- sum(muiPost>0)
	 whichNonZero_muiPost <- muiPost > 0
	 sum_muiPost <- sum(muiPost)
	 sumlog_muiPost <- sum(log(muiPost[whichNonZero_muiPost]))

	 nNonZero_muiAll <- nNonZero_muiPre +nNonZero_muiPost
	 sum_muiAll <- sum_muiPre+sum_muiPost
	 sumlog_muiAll <- sumlog_muiPre+sumlog_muiPost

	 # update overdispersion
	 phi <- update_phi(phi, mu, sum_muiAll, sumlog_muiAll, nNonZero_muiAll,v.phi=v.phi)

	 # update zero-inflation
	 psiB <- update_psiB(psiB,nNZ_b=nNonZero_muiPre, nZ_b=n.pre-nNonZero_muiPre, nNZ_a=nNonZero_muiPost,deltaPsi)
	 deltaPsi <- update_deltaPsi(deltaPsi, nNZ_b=0,nZ_b=0, nNZ_a=nNonZero_muiPost, psiB) 

	 # save samples
	 if(i>0 && i%%thin ==0){
	    idx <- i/thin
	    yPre.samples[idx,] <- yPre
	    yPost.samples[idx,] <- yPost
	    muiPre.samples[idx,] <- muiPre
	    muiPost.samples[idx,] <- muiPost
	    mu.samples[idx] <- mu
	    phi.samples[idx] <- phi
	    psiB.samples[idx] <- psiB
	    deltaPsi.samples[idx] <- deltaPsi
	 }

	 # save acceptance
	 accept <- accept + c(attributes(mu)$accept,attributes(phi)$accept, attributes(psiB)$accept, attributes(deltaPsi)$accept, attributes(muiPre)$accept, attributes(muiPost)$accept) 

	 # print status information
	 if(length(w <- which(niter.1percent %in% i))) {printMsg(w-1, time.start)}
      }
      # save last value of chain
      last <- list(yPre=yPre,yPost=yPost, mu=mu, muiPre=muiPre, muiPost=muiPost,phi=phi,psiB=psiB, deltaPsi=deltaPsi)
      # return results
      res <- list(mu=mu.samples, yPre=yPre.samples, yPost=yPost.samples, muiPre=muiPre.samples, muiPost=muiPost.samples, phi=phi.samples, psiPre=psiB.samples, deltaPsi=deltaPsi.samples, accept=accept/(nburnin+1+niter),last=last)
      return(res)
   }
   )
   #####


   result <- list(yPre=update_yPre, yPost=update_yPost,muiPre=update_muiPre,muiPost=update_muiPost, mu=update_mu, psiB=update_psiB, deltaPsi=update_deltaPsi, phi=update_phi, runMCMC=runMCMC) 
   return(result)
}

