debverbose <- function(string, level=1){
        if (.verboselevel >= level) cat(paste0( substr("           ",1, level+1),string,"\n"))
}


catchNArgamma <- function(n, shape, rate=1, minval=minrgamma){
    vector <- rgamma(n, shape, rate)
    if (sum(is.na(vector))>0)
        debverbose(paste0('Caught NA in rgamma with shape ',shape," and rate ",rate), level=5)
    vector[is.na(vector)] <- minval
    vector
}

divergenceKL <- function(a,b,c, invGamma=TRUE, logNormal=FALSE){
   #####################################
   KL.hg <- function(x, a,b,c, g="gamma"){
      logconstH <- log(2) +.5*(a-1)*{log(b)-log(c)} +log(besselK(2*sqrt(b*c),a-1))
      logdH <- (-a)*log(x) - c/x - b*x -logconstH
      dH <- exp(logdH)
      
      sq <- sqrt(a^2+4*b*c)
      
      logdG <- if(g=="gamma"){ 
        dgamma(x, shape=1+sq, rate=2*b +0.5*a*{a+sq}/c, log=TRUE)
      } else if(g=="invGamma"){
        dinvgamma(x, shape=sq-1, scale=2*c +0.5*a*{a-sq}/b,log=TRUE)
      } else {
	dlnorm(x,meanlog=1/sq +log(0.5*{-a+sq}/b),sdlog=(a^2+4*b*c)^(-.25),log=TRUE)
      }
      
      dH*(logdH-logdG)
   }
   ##########################################
   dKL  <- rep(Inf,3)
   d.hGa <- try(integrate(KL.hg, lower=0, upper=Inf,a=a,b=b,c=c,g="gamma")$value, silent=TRUE)
   if(!inherits(d.hGa,"try-error")) dKL[1]<- d.hGa

   if(invGamma){
      d.hIGa <- try(integrate(KL.hg, lower=0, upper=Inf,a=a,b=b,c=c, g="invGamma")$value,silent=TRUE)
      if(!inherits(d.hIGa,"try-error")) dKL[2] <- d.hIGa 
   }
   if(logNormal){ 
      d.hLN <- try(integrate(KL.hg, lower=0, upper=Inf,a=a,b=b,c=c, g="logNormal")$value,silent=TRUE)
      if(!inherits(d.hLN,"try-error")) dKL[3] <- d.hLN
   }
   # which is smaller?
#print(dKL)

   return(switch(which.min(dKL), "gamma","invGamma","logNormal")) # when in doubt use Gamma
   #return(dKL)
}

##' Update the parameter \eqn{theta} 
##' (with FC of form \eqn{h(x) = x^{-a} \exp(-bx -c/x)}{h(x) = x^(-a) exp(-bx -c/x)})) 
##' using eiter a gamma or an inverse gamma approximation to the FC as proposal
##' depending on which distribution has smaller KL divergence.
##'
##' @param argA argument a (in parts) of the function h (e.g. n*phi for theta=mu)
##' @param argC argument c of the function h (e.g. sum(mui)*phi for theta=mu)
##' @param a shape parameter of the gamma prior for \eqn{theta}
##' @param b rate parameter of the gamma prior for \eqn{theta}

##' @return updated parameter \eqn{theta}
##' @keywords internal
update_h_kl <- function(theta, argA, argC, a, b){ 
   # unnormalized log marginal posterior for theta

   
   log.post.theta <- function(theta){
      (-argA+a-1)*log(theta) - argC/theta - b*theta
   }

   #
   a1 <- argA-a+1
# cat("update_h_kl:",argA,argC,"\t") # todel
   if( is.na(a1)) a1 <- 0
   
   sq <- sqrt(a1^2+4*b*argC)

  
   ## decide which approximation to use
   # always use an  inverse gamma if argA-a+1 is > 2
   whichApprox <- "invGamma"
   # if argA-a+1 is <=1 then decide between a gamma and a logNormal approx
   #   when b>= 1 logNormal is mostly better
   if(a1 <= 1 ){
      if(b>=1){
	 whichApprox <- "logNormal"
      } else {
	 whichApprox <- divergenceKL(a=a1,b=b,c=argC, invGamma=FALSE,logNormal=TRUE)
      }
   } else if(a1 <=2){
      # for values of argA-a in (0,1] look at KL divergence (gamma, invGamma) to decide
      whichApprox <- divergenceKL(a=a1,b=b,c=argC, invGamma=TRUE,logNormal=TRUE)
   }

   theta.new <- switch(whichApprox,
      "invGamma"= MHstep(theta,rinvgamma, dinvgamma, log.post.theta, shape=sq-1, scale=2*argC +0.5*a1*{a1-sq}/b),
      "gamma"= MHstep(theta,rgamma, dgamma, log.post.theta, shape=1+sq, rate=2*b +0.5*a1*{a1+sq}/argC),
      "logNormal" = MHstep(theta,rlnorm, dlnorm, log.post.theta, meanlog=1/sq +log(0.5*{-a1+sq}/b),sdlog=(a1^2+4*b*argC)^(-.25))
   )

   debverbose(sprintf("Update theta(KL): %.3f->%.3f, using %s approximation",theta,theta.new,whichApprox), level=4)
   return(theta.new)
}


##' Update the parameter \eqn{theta} 
##' (with FC of form \eqn{h(x) = x^{-a} \exp(-bx -c/x)}{h(x) = x^(-a) exp(-bx -c/x)}) 
##' using a gamma approximation to the FC as proposal
##'
##' @param argA argument a (in parts) of the function h (e.g. n*phi for theta=mu)
##' @param argC argument c of the function h (e.g. sum(mui)*phi for theta=mu)
##' @param a shape parameter of the gamma prior for \eqn{theta}
##' @param b rate parameter of the gamma prior for \eqn{theta}

##' @return updated parameter \eqn{theta}
##' @keywords internal
update_h_gamma <- function(theta, argA, argC, a, b){  
   # unnormalized log marginal posterior for theta
   log.post.theta <- function(theta){
      (-argA+a-1)*log(theta) - argC/theta - b*theta
   }

   # compute mode & curvature
   mode.curv <- modeCurvature_h(a=argA-a+1, b=b, c=argC)
   # compute parameters of an approximating gamma distribution
   fc <- fc_approx_gamma(mode=mode.curv$mode,curvature=mode.curv$curvature)
   # propose new value
   theta.new <- MHstep(theta,rgamma, dgamma, log.post.theta, shape=fc["shape"], rate=fc["rate"])

   debverbose(sprintf("Update theta(gamma): %.3f->%.3f, using approximation shape %.3f and rate %.3f",
                      theta,theta.new, fc["shape"], rate=fc["rate"]), level=4)
   
   return(theta.new)
}

##' Update the parameter \eqn{theta} (with FC of form h) using an inverse gamma approximation to the FC as proposal
##'
##' @param argA argument a (in parts) of the function h (e.g. n*phi for theta=mu)
##' @param argC argument c of the function h (e.g. sum(mui)*phi for theta=mu)
##' @param a shape parameter of the gamma prior for \eqn{theta}
##' @param b rate parameter of the gamma prior for \eqn{theta}

##' @return updated parameter \eqn{theta}
##' @keywords internal
update_h_invGamma <- function(theta, argA, argC, a, b){ 
   # unnormalized log marginal posterior for theta
   log.post.theta <- function(theta){
      (-argA+a-1)*log(theta) - argC/theta - b*theta
   }

   # compute mode & curvature
   mode.curv <- modeCurvature_h(a=argA-a+1, b=b, c=argC)
#if(mode.curv$mode^2*mode.curv$curvature <1) cat("a =",argA-a+1,",\tb =",b,",\tc =",argC,"\n")
#if(argA <1.1) cat("mode=", mode.curv$mode,"; curv=",mode.curv$curvature,";\ta =",argA-a+1,";\tb =",b,";\tc =",argC,"\n")
   # compute parameters of an approximating gamma distribution
   fc <- fc_approx_invgamma(mode=mode.curv$mode,curvature=mode.curv$curvature)
   # propose new value
   theta.new <- MHstep(theta,rinvgamma, dinvgamma, log.post.theta, shape=fc["shape"], scale=fc["scale"])

   debverbose(sprintf("Update theta(invgamma): %.3f->%.3f, using approximation shape %.3f and rate %.3f",
                      theta,theta.new, fc["shape"], rate=fc["rate"]), level=4)

   
   return(theta.new)
}

##' Update the parameter \eqn{theta} (with FC of form h) using a lognormal approximation to the FC as proposal
##'
##' @param argA argument a (in parts) of the function h (e.g. n*phi for theta=mu)
##' @param argC argument c of the function h (e.g. sum(mui)*phi for theta=mu)
##' @param a shape parameter of the gamma prior for \eqn{theta}
##' @param b rate parameter of the gamma prior for \eqn{theta}

##' @return updated parameter \eqn{theta}
##' @keywords internal
update_h_logNormal <- function(theta, argA, argC, a, b){ 
   # unnormalized log marginal posterior for theta
   log.post.theta <- function(theta){
      (-argA+a-1)*log(theta) - argC/theta - b*theta
   }

   # compute mode & curvature
   mode.curv <- modeCurvature_h(a=argA-a+1, b=b, c=argC)
   # compute parameters of an approximating gamma distribution
   fc <- fc_approx_logNormal(mode=mode.curv$mode,curvature=mode.curv$curvature)
   # propose new value
   theta.new <- MHstep(theta,rlnorm, dlnorm, log.post.theta, meanlog=fc["mean"], sdlog=fc["sd"])

   debverbose(sprintf("Update theta(logNormal): %.3f->%.3f, using approximation meanlog %.3f and sdlog %.3f",
                      theta,theta.new, fc["mean"], rate=fc["sd"]), level=4)

   
   return(theta.new)
}


# ---------------------
##' Update latent egg counts \eqn{y_1,\ldots,y_n} based on a truncated Poisson distribution
##'
##' @param mui mean
##' @param n number of samples
##' @param propC 1 minus 1/f where f is the dilution factor
##' @param fec observed number of eggs (on the McMaster slide)

##' @return updated parameter vector \eqn{y}
##' @keywords internal
update_Y <- function(mui, n, propC, fec){
   rpois(n, lambda = propC*mui) + fec
}

#-------------------------
##' Update the overdispersion parameter \eqn{phi} 
##' using a uniform random walk proposal
##'
##' @param phi phi
##' @param mu mu
##' @param sum_mui sum of \eqn{mu_i}
##' @param sumlog_mui sum of log(\eqn{mu_i})
##' @param a shape parameter of the gamma prior for \eqn{phi}
##' @param b rate parameter of the gamma prior for \eqn{phi}

##' @return updated parameter \eqn{phi}
##' @keywords internal
update_phi_unif <- function(phi,mu, sum_mui,sumlog_mui, n, a, b, v.phi){  
   # unnormalized log marginal posterior for phi
   log.post.phi <- function(phi){
      (n*phi+a-1)*log(phi) - n*lgamma(phi) + phi*(-n*log(mu)+ sumlog_mui) - sum_mui*phi/mu - b*phi  
   }

   # do an MH step
   phi.new <- MH_RW_unif(phi,log.post.phi, v=v.phi)

   debverbose(sprintf("Update phi(unif): %.3f->%.3f",phi,phi.new), level=3)

   
   return(phi.new)
}

#~~~~~~~~~~~~~~~~~~~~~++
update_phiTransf_unif <- function(phi,mu, sum_mui,sumlog_mui, n, a, b, v.phi, factor=5){  

#cat("phi=",phi,"phi01=",atan(phi)*2/pi," mu=",mu," sumMui=",sum_mui,"\n")
   phi.01 <- atan(phi/factor)*2/pi

   # unnormalized log marginal posterior for phi
   log.post.phi <- function(phi.01){
      phi <- factor*tan(phi.01*pi/2)
      (n*phi)*log(phi) - n*lgamma(phi) + phi*(-n*log(mu)+ sumlog_mui) - sum_mui*phi/mu +(a-1)+log(phi.01) + (b-1)*log(1-phi.01) 
   }

   # do an MH step
   phi.new <- MH_RW_unif01(phi.01,log.post.phi, v=v.phi) # truncate at 0 and 1? use a beta as proposal?
   
   phi.new <- factor*tan(phi.new*pi/2)

   debverbose(sprintf("Update phi.transf(unif): %.3f->%.3f",phi.01,phi.new), level=3)

   return(phi.new)
}

 
##' Update the overdispersion parameter \eqn{phi} 
##' using a gamma approximation to the FC as proposal
##'
##' @param phi phi
##' @param mu mu
##' @param sum_mui sum of \eqn{mu_i}
##' @param sumlog_mui sum of log(\eqn{mu_i})
##' @param a shape parameter of the gamma prior for \eqn{phi}
##' @param b rate parameter of the gamma prior for \eqn{phi}

##' @return updated parameter \eqn{phi}
##' @keywords internal
update_phi_gamma <- function(phi,mu, sum_mui,sumlog_mui, n, a, b,...){  
   # unnormalized log marginal posterior for phi
   log.post.phi <- function(phi){
      (n*phi+a-1)*log(phi) - n*lgamma(phi) + phi*(-n*log(mu)+ sumlog_mui) - sum_mui*phi/mu - b*phi  
   }

   # compute mode & curvature
   mode.curv <- modeCurvature_phi(n,mu,logprodmu=sumlog_mui, summu.mu=sum_mui/mu,a=a, b=b)
   # compute parameters of an approximating gamma distribution
   fc <- fc_approx_gamma(mode=mode.curv$mode,curvature=mode.curv$curvature)
   # do an MH step
   phi.new <- MHstep(phi,rgamma,dgamma, log.post.phi, shape=fc["shape"], rate=fc["rate"])

   debverbose(sprintf("Update phi(gamma): %.3f->%.3f",phi,phi.new), level=3)
   
   return(phi.new)
}
#-------------------------
##' Update a parameter vector theta \eqn{theta_1,\ldots,theta_n} based on a gamma distribution
##' using a gamma proposal with slightly shifted shape parameter to avoid to small theta
##'
##' @param theta theta
##' @param s shape
##' @param r rate
##' @param n length of theta

##' @return updated parameter \eqn{theta}
##' @keywords internal
update_theta_gammaAdd <- function(theta, s,r, n){
   addS <- 0.05
   # draw from full conditional if shape is not "too small"
   if(all(s>0.05)){
#      ret <- catchNA( rgamma(n, shape = s, rate= r), minrgamma)
      ret <- catchNArgamma(n, shape = s, rate= r, minrgamma)
      attr(ret,"accept") <- rep(TRUE, n)
   } else { # shift scale a bit to avoid too small theta's
      # propose a new value
      addS <- ifelse(s<=0.05,addS,0)
#      new <- catchNA( rgamma(n,shape=s+addS, rate=r), minrgamma)
      new <- catchNArgamma(n,shape=s+addS, rate=r, minrgamma)
      # compute difference of log(posterior ratio)
      diff.log.post <- dgamma(new, shape=s, rate=r, log=TRUE) - dgamma(theta, shape=s, rate=r, log=TRUE)
      # compute difference of log(proposal ratio)
      diff.log.prop <- dgamma(theta, shape=s+addS, rate=r, log=TRUE) - dgamma(new, shape=s+addS, rate=r, log=TRUE)
      # compute acceptance ratio
      acc.ratio <- exp(min(0, diff.log.post + diff.log.prop))

      u <- runif(n)
      accept <- u < acc.ratio
      ret <- ifelse(accept, new, theta)
      
      attr(ret,"accept") <- accept      
   }
   return(ret)
}
#-------------------------
##' Update a parameter vector theta \eqn{theta_1,\ldots,theta_n} based on a gamma distribution
##'
##' @param s shape
##' @param r rate
##' @param n length of theta

##' @return updated parameter \eqn{theta}
##' @keywords internal
update_theta_gamma <- function(s,r, n){
   catchNArgamma(n, shape = s, rate= r, minrgamma)
#      catchNA( rgamma(n, shape = s, rate= r), minrgamma)
}

##' Update a parameter \eqn{theta} based on a gamma distribution
##'
##' @param s s
##' @param r r
##' @param a shape parameter of the gamma prior for \eqn{theta}
##' @param b rate parameter of the gamma prior for \eqn{theta}
##' @param \dots extra arguments (not used)

##' @return updated parameter \eqn{theta}
##' @keywords internal
update_theta_gamma1 <- function(s,r, a,b, ...){
#   catchNA( rgamma(1, shape = s+a, rate= r+b), minrgamma)
   catchNArgamma(1, shape = s+a, rate= r+b, minrgamma)
}


##' Update a parameter vector theta \eqn{theta_1,\ldots,theta_n} based on a gamma mixture distribution
##'
##' @param mui mui
##' @param y y
##' @param mu mu
##' @param phi phi
##' @param psi mixture probability
##' @param n length of theta
##' @param nNonZero_fec number of non-zero observed counts
##' @param whichZero_fec vector of logicals indicating which observed counts are zero
##' @param p p, defaults to 1

##' @return updated parameter \eqn{theta}
##' @keywords internal
update_theta_gammaMix <- function(mui, y, mu, phi, psi, n, whichZero_fec, nNonZero_fec, p=1){

#cat("upd_theta_gammaMix: mui=",mui," y=", y, " mu=",mu, " phi=",phi, " psi=",psi, "\n")

   # propose new values
   # in the case of non-zero observations, this corresponds to a draw from the full conditional
   # in the case of a zero observation, an MH step is needed
   # the proposal distribution is a mixture (1-psi)*I(mui=0) + psi*Ga(phi, 1+phi/mu)
   # if phi is very small, this leads to proposed mui's approx=0
   # -> add 0.1 to the shape of proposal in case of zero observations to avoid very small proposed mui's
   addS <- 0.1
   new <- numeric(n)

   # draw mui's for non-zero observed counts from a gamma distribution
   shape <- y+phi
   rate <- p+phi/mu
   # for non-zero observed counts, draw from gamma full conditional
#   new[!whichZero_fec] <- catchNA( rgamma(nNonZero_fec, shape=shape[!whichZero_fec], rate=rate), minrgamma)
   new[!whichZero_fec] <- catchNArgamma(nNonZero_fec, shape=shape[!whichZero_fec], rate=rate, minrgamma)

if(n==nNonZero_fec) return(new)

   yZero_fec <- y[whichZero_fec]
   muiZero_fec <- mui[whichZero_fec]
   #propose mui values for zero observed counts
   proposed_mui <- new[whichZero_fec]

   # draw mui's for zero observed counts from a mixture of a gamma distribution and a point mass at zero
   which_zeroFEC_posMui <- rbinom(n-nNonZero_fec, size=1, prob=psi)==1
#   proposed_mui[which_zeroFEC_posMui] <- catchNA(
#       rgamma(sum(which_zeroFEC_posMui),shape=yZero_fec[which_zeroFEC_posMui]+phi+addS, rate=rate), minrgamma)
   proposed_mui[which_zeroFEC_posMui] <- catchNArgamma(sum(which_zeroFEC_posMui),shape=yZero_fec[which_zeroFEC_posMui]+phi+addS, rate=rate, minrgamma)


   nZero <- sum(whichZero_fec)

   log.post.muiZero <- function(mui){
      ifelse(mui==0,
		  log((1-psi)*(yZero_fec==0)),
		  log(psi) + dgamma(mui,shape=phi, rate=phi/mu, log=TRUE) + dpois(yZero_fec, lambda=mui,log=TRUE)
	    )
	    
   }
   
   # compute difference of log(proposal ratio) for zero fec counts
   diff.log.prop <- logdgammaMixture(muiZero_fec, shape=yZero_fec+phi+addS, rate=rate, prob=psi)- logdgammaMixture(proposed_mui, shape=yZero_fec+phi+addS, rate=rate, prob=psi)
   # compute difference of log(posterior ratio) for zero obserations
   diff.log.post <- log.post.muiZero(proposed_mui) -log.post.muiZero(muiZero_fec)
   # accept value? 
   acc.ratio <- exp(pmin(0,diff.log.post+diff.log.prop))
   u <- runif(nZero)
   accept <- u < acc.ratio

#cat("y: nZero",nZero,"new=",new,"AR=",acc.ratio,"u",u,"accept",accept,"newZero",newZero,"yZero",yZero,"\n")
#cat("y: nZero",nZero,"yZero",yZero,"newZero",newZero,"AR=",acc.ratio,"u",u,"\n")

   new[whichZero_fec] <- ifelse(accept, proposed_mui, muiZero_fec)
   
   attr(new,"accept") <- accept
   return(new)

}



#----------------------------------
##' Update a parameter \eqn{theta} based on a beta distribution
##'
##' @param x number of events
##' @param n sample size
##' @param a shape1 parameter of the beta prior for \eqn{theta}
##' @param b shape2 parameter of the beta prior for \eqn{theta}

##' @return updated parameter \eqn{theta}
##' @keywords internal
update_theta_beta <- function(x,n, a,b){
#cat("upd_theta_beta: shape1=",x+a," shape2=", n-x+b, "\n")
   rbeta(1, shape1 = x+a, shape2= n-x+b)
}

#----------------------------------------
##' Update the zero-inflation parameter \eqn{psi} 
##' using a beta approximation to the FC as proposal
##'
##' @param psi psi
##' @param nNZ_b nNZ.b
##' @param nZ_b nZ.b
##' @param nNZ_a nNZ.a
##' @param delta delta
##' @param n n
##' @param a shape1 parameter of the beta prior for \eqn{psi}
##' @param b shape2 parameter of the beta prior for \eqn{psi}

##' @return updated parameter \eqn{phi}
##' @keywords internal
update_psi_beta <- function(psi,nNZ_b,nZ_b,nNZ_a,delta, n, a, b ){
   # unnormalized log marginal posterior for psiB
   log.post.psi <- function(psi){
      (n-nNZ_a)*log(1-delta*psi) + (nZ_b+b-1)*log(1-psi) + (nNZ_b+nNZ_a +a-1)*log(psi)
   }
   
   # compute mode & curvature
   mode.curv <- modeCurvature_psi(d=delta, a=n-nNZ_a, b=nZ_b+b-1, c=nNZ_b+nNZ_a+a-1)
   # compute parameters of an approximating beta distribution
   fc <- fc_approx_beta(mode=mode.curv$mode,curvature=mode.curv$curvature, dens01=mode.curv$dens01)
   # propose new value
   psi.new <- MHstep(psi,rbeta, dbeta, log.post.psi, shape1=fc["shape1"], shape2=fc["shape2"])


   
   return(psi.new)
}

#----------------------------------------
##' Update the reduction in mean parameter \eqn{delta} 
##' using a beta approximation to the FC as proposal
##'
##' @param delta delta
##' @param sumYpost sumYpost
##' @param sumMuPre sumMuPre
##' @param a shape1 parameter of the beta prior for \eqn{delta}
##' @param b shape2 parameter of the beta prior for \eqn{delta}

##' @return updated parameter \eqn{delta}
##' @keywords internal
update_delta_beta <- function(delta, sumYpost, sumMuPre, a, b ){
   # unnormalized log marginal posterior for delta
   log.post.delta <- function(delta){
      (sumYpost+a-1)*log(delta) + (b-1)*log(1-delta) - delta*sumMuPre
   }
   
   # compute mode & curvature
   mode.curv <- modeCurvature_delta(a=sumYpost+a-1, b=b-1, c=sumMuPre)
   # compute parameters of an approximating beta distribution
   fc <- fc_approx_beta(mode=mode.curv$mode,curvature=mode.curv$curvature, dens01=NULL)
   # propose new value
   delta.new <- MHstep(delta,rbeta, dbeta, log.post.delta, shape1=fc["shape1"], shape2=fc["shape2"])
   
   return(delta.new)
}

