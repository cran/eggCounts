##' Convert the list with samples to a \code{mcmc} object
##'
##' @param x result of a call to \code{fec_mcmc} or \code{fecr_mcmc}
##' @return list with \code{mcmc} objects

##' @keywords internal
##' @export
samples2mcmc <- function(x){
   model <- x$model
   samples <- x$samples
   start <- x$nburnin + 1
   thin <- x$thin

   #################################
   ## one-sample situation: (ZI)PoGa
   ################################
   if(model%in% c("PoGa","ZIPoGa","PoGa2","ZIPoGa2")){
      mu <- samples$mu
      phi <- samples$phi
      psi <- samples$psi # NULL in case of ZIPoGa

      mui <- samples$mui
      colnames(mui) <- paste("mui",1:ncol(mui),sep="")

      if(model %in% c("PoGa2","ZIPoGa2")){
	 y <- NULL
      } else {
      y <- samples$y
      colnames(y) <- paste("y",1:ncol(y),sep="")
      y <- mcmc(y, start=start, thin=thin)
      }

      all1 <- mcmc(cbind(mu=mu, phi=phi,psi=psi),start=start, thin=thin)
      FEC <- mcmc(cbind(meanEpG=rowMeans(mui)),start=start, thin=thin)
      mui <- mcmc(mui, start=start, thin=thin)
      
      return(list(fec=FEC,all1=all1,mui=mui, y=y))

   } else {
   #################################
   ## two-sample situation
   #################################
      mu <- samples$mu
      phi <- samples$phi
      
      muiPre <- samples$muiPre
      colnames(muiPre) <- paste("muiPre",1:ncol(muiPre),sep="")
      muiBar.pre <-rowMeans(muiPre)
      muiC <- mcmc(muiPre, start=start, thin=thin)

      yPre <- samples$yPre
      colnames(yPre) <- paste("yPre",1:ncol(yPre),sep="")
      yC <- mcmc(yPre, start=start, thin=thin)

      yPost <- samples$yPost
      colnames(yPost) <- paste("yPre",1:ncol(yPost),sep="")
      yT <- mcmc(yPost, start=start, thin=thin)


      if(model=="paired"){
	 deltaMu <- samples$delta
	 # compute variables of interest
	 muiPost <- muiPre*deltaMu  # matrix * vector
         colnames(muiPost) <- paste("muiPost",1:ncol(muiPost),sep="")
	 muiBar.post <- rowMeans(muiPost)
	 fecr <- (1-deltaMu)

	 # all parameters in the model (except mu_i, y_i) + computed params
	 all1 <- mcmc(cbind(mu=mu, phi=phi, deltaMu=deltaMu), start=start, thin=thin)
	 FECR <- mcmc(cbind(fecr=fecr, meanEpg.untreated=muiBar.pre, meanEpg.treated=muiBar.post), start=start, thin=thin)
	 muiT <- mcmc(muiPost, start=start, thin=thin)
	 
      } else if(model=="unpaired"){
	 deltaMu <- samples$delta

	 muiPost <- samples$muiPost
	 colnames(muiPost) <- paste("muiPost",1:ncol(muiPost),sep="")
	 muiBar.post <- rowMeans(muiPost)
	 fecr <- 1-muiBar.post/muiBar.pre

	 all1 <- mcmc(cbind(mu=mu, phi=phi,deltaMu=deltaMu), start=start, thin=thin)
	 FECR <- mcmc(cbind(fecr=fecr, meanEpg.untreated=muiBar.pre, meanEpg.treated=muiBar.post), start=start, thin=thin)
	 muiT <- mcmc(muiPost, start=start, thin=thin)
	 
      } else if (model=="ZIPoGa_u"){
	 psiPre <- samples$psiPre
	 psiPost <- samples$psiPost

	 muiPost <- samples$muiPost
	 colnames(muiPost) <- paste("muiPost",1:ncol(muiPost),sep="")
	 muiBar.post <- rowMeans(muiPost)
 	 fecr <- 1-muiBar.post/muiBar.pre
     
	 all1 <- mcmc(cbind(mu=mu, phi=phi,psiTreated=psiPre, psiControl=psiPost), start=start, thin=thin)
	 FECR <- mcmc(cbind(fecr=fecr, meanEpg.untreated=muiBar.pre, meanEpg.treated=muiBar.post), start=start, thin=thin)
	 muiT <- mcmc(muiPost, start=start, thin=thin)
	 
      } else if (model=="ZIPoGa_upd"){
	 psiPre <- samples$psiPre
	 deltaPsi <- samples$deltaPsi
	 psiPost <- psiPre*deltaPsi

	 muiPost <- samples$muiPost
	 colnames(muiPost) <- paste("muiPost",1:ncol(muiPost),sep="")
	 muiBar.post <- rowMeans(muiPost)
 	 fecr <- 1-muiBar.post/muiBar.pre
     
	 all1 <- mcmc(cbind(mu=mu, phi=phi, deltaPsi=deltaPsi,psiTreated=psiPre, psiControl=psiPost), start=start, thin=thin)
	 FECR <- mcmc(cbind(fecr=fecr, meanEpg.untreated=muiBar.pre, meanEpg.treated=muiBar.post), start=start, thin=thin)
	 muiT <- mcmc(muiPost, start=start, thin=thin)
	 
      }

      if ((model!="paired")&(min(muiBar.pre) <= .Machine$double.eps))
          cat("NOTE: zero before treatment means appear.\n      You need to expect additional warnings/errors.\n") 
      
      return(list(fecr=FECR, all1=all1,muiControl=muiC,muiTreated=muiT, yControl=yC,yTreated=yT))

   }

}

##' Print information about mcmc run
##'
##' @param x result of a call to \code{fec_mcmc}
##' @param digits digits passed to printing function
##' @param quantiles a vector of quantiles to evaluate for each variable
##' @return x

##' @keywords internal
##' @method print fecm
##' @author Michaela Paul, with contributions from Reinhard Furrer 
##' @export
print.fecm <- function(x, digits=3, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), ...){
   samples <- samples2mcmc(x)
   cat("Model: ", x$model,"\n")
   #cat("Number of samples: ",x$nsamples," (burnin=",x$nburnin,", thinning=",x$thin,")\n", sep="")
   print(summary(samples$all, quantiles=quantiles), digits=digits)
   invisible(x)
}

##' Print information about mcmc run
##'
##' @param x result of a call to \code{fecr_mcmc}
##' @param digits digits passed to printing function
##' @param quantiles a vector of quantiles to evaluate for each variable
##' @return x

##' @keywords internal
##' @author Michaela Paul, with contributions from Reinhard Furrer 
##' @method print fecrm
##' @export
print.fecrm <- function(x, digits=3, quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975), ...){
   samples <- samples2mcmc(x)
   cat("Model: ", x$model,"\n")
   #cat("Number of samples: ",x$nsamples," (burnin=",x$nburnin,", thinning=",x$thin,")\n", sep="")
   print(summary(samples$fecr, quantiles=quantiles), digits=digits)
   invisible(x)
}
