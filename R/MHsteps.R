checkMHiterpars <- function( maxiter.pilot,  nburnin, nsamples, thin)
{
    if ((maxiter.pilot<0)|(ceiling(maxiter.pilot)!=floor(maxiter.pilot)))
        stop("'maxiter.pilot' should be a positive integer", call.=FALSE)
    if ((nburnin < 0)|(ceiling(nburnin)!=floor(nburnin)))
        stop("'nburnin' should be a positive integer", call.=FALSE)
    if ((nsamples < 0)|(ceiling(nsamples)!=floor(nsamples)))
        stop("'nburnin' should be a positive integer", call.=FALSE)
    if ((thin < 1)|(ceiling(thin)!=floor(thin)))
        stop("'thin' should be strictly a positive integer", call.=FALSE)    
    
    if (maxiter.pilot < 2)  cat("NOTE: 'maxiter.pilot' seems low.\n")
    if (nburnin < 100)      cat("NOTE: 'nburnin' seems low.\n")
    if (nsamples < 100)     cat("NOTE: 'nsamples' seems low.\n")
    
    invisible()
}
   
    
    

##' Check acceptance rates 
##' 
##' The function determines whether the acceptance rate of a given vector with 
##' acceptance states lies between 0.35 and 0.45
##' 
##' @param accept 0/1 vector indicating acceptance
##' @param v used tuning parameter
##' @param verbose print progress information

##' @return vector with value of tuning parameter and a boolean indicating an acceptable rate
##'
##' @keywords internal
checkAR <- function(accept, v, verbose=TRUE){  ## todo accept==AR
   AR.ok <- FALSE
   # compute acceptance rate
   AR <- mean(accept)
   if(AR < 0.29){
      v.new <- v*0.75
      if(verbose) cat(" decrease tuning parameter", sQuote(substitute(v)), "to",round(v.new,3),"\n")
   } else if (AR >0.40){
      v.new <- v*1.25
      if(verbose) cat(" increase tuning parameter", sQuote(substitute(v)), "to",round(v.new,3),"\n")
   } else {
      AR.ok <- TRUE
      v.new <- v
      #cat("tuning parameter", sQuote(substitute(v)), "selected as:",v,"\n")
   }	 
   return(c(v.new,AR.ok))
}


##' Do a Metropolis-Hastings step
##' 
##' @param old current state of the chain
##' @param rprop a random number generating function with arguments n=number of random variables, ... extra-arguments
##' @param dprop the corresponding density function with argument x, and the same ... arguments
##' @param logUpost function evaluating the log unnormalized marginal posterior
##' @param \dots extra arguments for the proposal density and random number generator

##' @return new state of the chain with a boolean indicating acceptance as attribute
##'
##' @keywords internal
MHstep <- function(old,rprop,dprop, logUpost,...){
   # propose a new value
   new <- rprop(1,...)
   if(!is.finite(new)){
      ret <- old
      attr(ret,"accept") <- FALSE
      return(ret)
   }
   # compute difference of log(posterior ratio)
   diff.log.post <- logUpost(new) - logUpost(old)
   # compute difference of log(proposal ratio)
   diff.log.prop <- dprop(x=old,...,log=TRUE) - dprop(x=new,...,log=TRUE)
   # compute acceptance ratio
   acc.ratio <- exp(min(0, diff.log.post + diff.log.prop, na.rm=TRUE))

   u <- runif(1)
   accept <- u < acc.ratio
   ret <- ifelse(accept, new, old)
   
   attr(ret,"accept") <- accept

   debverbose(sprintf("MH step: %.3f->%.3f (%d); diff post: %.3f, diff prop: %.3f, AR: %.3f",
                      old, new, accept,   diff.log.post,diff.log.prop,acc.ratio), level=3)
   
   return(ret)
}

##' MH step with uniform (truncated at 0) proposal around the current value
##' 
##' @param old current state of the chain
##' @param logUpost function evaluating the log unnormalized marginal posterior
##' @param v tuning parameter, detemines the width of the uniform proposal

##' @return new state of the chain with a boolean indicating acceptance as attribute
##'
##' @keywords internal
MH_RW_unif <- function(old, logUpost, v){
   # propose new value
   new <- runif(1, min = max(old-v, 0), max = old +v)

   # compute difference of log(posterior ratio)
   diff.log.post <- logUpost(new) - logUpost(old)

   # compute difference of log(proposal ratio)
   diff.log.prop <- -log(old + v - max(old - v,0)) + log(new + v -max(new - v, 0))

   acc.ratio <- exp(min(0, diff.log.post + diff.log.prop, na.rm=TRUE))

   u <- runif(1)
   accept <- u < acc.ratio
   ret <- ifelse(accept, new, old)
   
   attr(ret,"accept") <- accept

   debverbose(sprintf("MH step(unif): %.3f->%.3f (%d); diff post: %.3f, diff prop: %.3f, AR: %.3f",
                      old, new, accept,   diff.log.post,diff.log.prop,acc.ratio), level=3)

   
   return(ret)

}

##' MH step with uniform (truncated at 0 and 1) proposal around the current value
##' 
##' @param old current state of the chain
##' @param logUpost function evaluating the log unnormalized marginal posterior
##' @param v tuning parameter, detemines the width of the uniform proposal

##' @return new state of the chain with a boolean indicating acceptance as attribute
##'
##' @keywords internal
MH_RW_unif01 <- function(old, logUpost, v){
   # propose new value
   new <- runif(1, min = max(old-v, 0), max = min(old +v,1))

   # compute difference of log(posterior ratio)
   diff.log.post <- logUpost(new) - logUpost(old)

   # compute difference of log(proposal ratio)
   diff.log.prop <- -log(min(old + v,1) - max(old - v,0)) + log(min(new + v,1) -max(new - v, 0))

   acc.ratio <- exp(min(0, diff.log.post + diff.log.prop, na.rn=TRUE))

   u <- runif(1)
   accept <- u < acc.ratio
   ret <- ifelse(accept, new, old)
   
   attr(ret,"accept") <- accept
   debverbose(sprintf("MH step(unif01): %.3f->%.3f (%d); diff post: %.3f, diff prop: %.3f, AR: %.3f",
                      old, new, accept,   diff.log.post,diff.log.prop,acc.ratio), level=3)

   return(ret)

}

##' Approximate a full conditional with given mode and curvature by a gamma distribution
##'
##' @param mode mode of the full conditional
##' @param curvature curvature at the mode
##' @return shape and rate parameter of the approximating gamma distribution
##'
##' @keywords internal
fc_approx_gamma <- function(mode, curvature){
   # rate parameter
   beta <- curvature*mode
   # shape parameter
   alpha <- mode*beta +1
   return(c(shape=alpha,rate=beta))   
}


##' Approximate a full conditional with given mode and curvature by an inverse gamma distribution
##'
##' @param mode mode of the full conditional
##' @param curvature curvature at the mode
##' @return shape and scale parameter of the approximating inverse gamma distribution
##'
##' @keywords internal
fc_approx_invgamma <- function(mode, curvature){
   # scale parameter
   beta <- curvature*mode^3
   # shape parameter
   alpha <- beta/mode-1

   if(alpha <=0){
    cat("NOTE: shape parameter of InvGamma approximation negative!,",alpha," -> set to 1\n")
#cat("mode=",mode, "curvature=",curvature,"curv*mode^2 =",curvature*mode^2,"\n")
      alpha <- 1
   }
   return(c(shape=alpha,scale=beta))   
}

##' Approximate a full conditional with given mode and curvature by a log normal distribution
##'
##' @param mode mode of the full conditional
##' @param curvature curvature at the mode
##' @return mean and standard deviation (on log-scale) of the approximating log normal distribution
##'
##' @keywords internal
fc_approx_logNormal <- function(mode, curvature){
   # mean
   mu <- log(mode) +1/(curvature*mode^2)
   # variance
   sigma2 <- 1/(curvature*mode^2)
   return(c(mean=mu,sd=sqrt(sigma2)))
}

##' Approximate a full conditional with given mode and curvature by a beta distribution
##'
##' @param mode mode of the full conditional
##' @param curvature curvature at the mode
##' @param dens01 density at the boundary
##' @return shape1 and shape2 parameter of the approximating beta distribution
##'
##' @keywords internal
fc_approx_beta <- function(mode, curvature, dens01){
   if(!is.null(curvature)){
      alpha <-  -curvature*mode^3+curvature*mode^2+1
      beta <- curvature*mode^3-2*curvature*mode^2+curvature*mode+1
   } else if(round(mode,10)==1){
      beta <- 1
      alpha <- dens01
   } else if(round(mode,10)==0){
      alpha <- 1
      beta <- dens01
   }
   return(c(shape1=alpha,shape2=beta))
}


##' Compute mode and curvature at the mode of a full conditional of form
##' \eqn{h(x) = x^{-a} \exp(-bx -c/x)}{h(x) = x^(-a) exp(-bx -c/x)}
##'
##' @param a a
##' @param b b
##' @param c c

##' @return vector with mode and curvature at the mode
##'
##' @keywords internal
modeCurvature_h <- function(a,b,c){
   mode <- 0.5*(-a+sqrt(a^2+4*b*c))/b
   curvature <- -a/mode^2+2*c/mode^3
   return(list(mode=mode, curvature=curvature))  
}


##' Compute mode and curvature at the mode of a full conditional of form
##' \eqn{h(x) = (1-d x)^{a} (1-x)^{b} x^{c}}{h(x) = (1-d*x)^a (1-x)^b x^c}
##'
##' @param a a
##' @param b b
##' @param c c
##' @param d d

##' @return list with mode, curvature at the mode, and the density const*h(mode) should the mode be at the boundary
##'
##' @keywords internal
modeCurvature_psi <- function(a,b,c,d){
   m <-  (b + c + a*d + c*d - sqrt(b^2 + 2*b*(c + a*d - c*d) + (c - (a + c)*d)^2))/(2*(a + b + c)*d)
   
   # if mode is 0 or 1, then compute density at point 0 (or 1)
   curvature <- dens01 <- NULL
   if(round(m,10) ==1 | round(m,10)==0){
      # - log(h)
      h_psi <- function(x, a,b,c,d){ (1-d*x)^a*(1-x)^b*x^c }
      dens01 <- h_psi(m, a=a, b=b,c=c, d=d) / integrate(h_psi,lower=0,upper=1,d=d, a=a,b=b,c=c)$value
   } else {
      # 2nd deriv at mode
      curvature <- a*d^2/(1-d*m)^2+b/(1-m)^2+c/m^2
   }
  return(list(mode=m,curvature=curvature,dens01=dens01))
}

##' Compute mode and curvature at the mode of a full conditional of form
##' \eqn{h(x) = x^{a} (1-x)^{b} \exp(-c x)}{h(x) = x^a (1-x)^b exp(c*x)}
##'
##' @param a a
##' @param b b
##' @param c c

##' @return list with mode, curvature at the mode, and the density const*h(mode) should the mode be at the boundary
##'
##' @keywords internal
modeCurvature_delta <- function(a,b,c){
   m <- -(1/2)*(-a-b-c+sqrt(a^2+2*a*b-2*c*a+b^2+2*b*c+c^2))/c
  
   curvature <- a/m^2+b/(1-m)^2

  return(list(mode=m,curvature=curvature))
}

##' Numerically compute mode and curvature at the mode of a full conditional for \eqn{\phi}{phi} of form
##' -(n*phi+a-1)*log(phi) + n*lgamma(phi) + n*phi*log(mu)- (phi-1)*logprodmu + phi*(summu.mu+b)
##'
##' @param n n
##' @param mu mu
##' @param logprodmu logprodmu
##' @param summu.mu summu.mu
##' @param a a
##' @param b b

##' @return list with mode and curvature at the mode
##'
##' @keywords internal
modeCurvature_phi <- function(n,mu,logprodmu,summu.mu,a,b){
   # - log full conditional for phi
   logfc <- function(phi) -(n*phi+a-1)*log(phi) + n*lgamma(phi) + n*phi*log(mu)- phi*logprodmu + phi*(summu.mu+b)  
   # first derivative
   logfc.d1 <- function(phi) -n*log(phi)-(n*phi+a-1)/phi+n*digamma(phi)+n*log(mu)-logprodmu+summu.mu +b

   # find mode using uniroot (seems to be slightly faster than optimize)
   # phi is usually in the range between  0 and 10
   mode <- uniroot(logfc.d1,lower=1e-10,upper=100,tol=1e-3)$root
#    m <- uniroot(logfc.d1,lower=1e-10,upper=100,tol=1e-3)
# cat("phi uniroot iterations",m$iter,"\n")
#    mode <- m$root
   #mode <- optimize(logfc, lower=1e-10,upper=1000)
   curvature <- -2*n/mode+(n*mode+a-1)/mode^2+n*trigamma(mode)

   return(list(mode=mode, curvature=curvature))   
}

