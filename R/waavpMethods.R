######################################
## compute standard FECR
#######################################

##' Compute standard FECRT according to WAAVP guidelines
##'
##' Computes the standard Faecal Egg Count Reduction test together with approximate
##' confidence intervals according to the WAAVP guidelines 
##' (Coles et al., 1992, 2006). The function also returns bootstrap percentile 
##' confidence intervals.

##' @param epg1 faecal egg counts in untreated animals
##' @param epg2 faecal egg counts in treated animals
##' @param paired logical indication whether or not the samples are paired 
##' @param alpha confidence level of the intervals
##' @param R number of bootstrap replicates
##' @param \dots extra arguments (not used)

##' @return A list with
##' \item{estimate}{the estimated percentage reduction in mean epg rate}
##' \item{bootCI}{corresponding percentile bootstrap confidence interval}
##' \item{approxCI}{corresponding approximate confidence interval}

##' @references 
##' Coles GC, Bauer C, Borgsteede FHM, Geerts S, Klei TR, Taylor MA,  Waller, PJ (1992).
##' World Association for the Advancement of Veterinary Parasitology (WAAVP)
##' methods for the detection of anthelmintic resistance in nematodes of 
##' veterinary importance, Veterinary Parasitology, 44:35-44.
##'
##' Coles GC, Jackson F, Pomroy WE, Prichard RK,  von Samson-Himmelstjerna G,
##' Silvestre A, Taylor MA,  Vercruysse J (2006).
##' The detection of anthelmintic resistance in nematodes of veterinary	importance,
##' VeterinaryParasitology, 136:167-185.

##' @examples
##' data(epgs)
##' fecrtCI(epgs$before, epgs$after, paired=TRUE)

##' @export
# Coles:       FECRT = 1 - T2/C2   -> paired = FALSE
# Kochapakdee: FECRT = 1 - T2/T1   -> paired = TRUE
fecrtCI <- function(epg1, epg2, paired=FALSE, alpha=0.05, R=1999,...){

   conf.level <- 1-alpha

   if(paired){
      nas <- is.na(epg1) | is.na(epg2)
      epg1 <- epg1[!nas]
      epg2 <- epg2[!nas]
   } else {
      epg1 <- epg1[!is.na(epg1)]
      epg2 <- epg2[!is.na(epg2)]
   }

   ni <- c(length(epg1), length(epg2))
   data <- data.frame(epgs = c(epg1,epg2), grp = as.factor(rep(1:2, ni)))
   meanratio_unpaired <- function(d, i) {
      ind1 <- i[1:ni[1]]
      ind2 <- i[-(1:ni[1])]
      x1 <- d[ind1, 1]
      x2 <- d[ind2, 1]
      1-mean(x2)/mean(x1)
   }
   meanratio_paired <- function(d,i){
      ind1 <- i[1:ni[1]]
      ind2 <- ind1+ni[1]
      x1 <- d[ind1, 1]
      x2 <- d[ind2, 1]
      1-mean(x2)/mean(x1)
   }

   if(paired){
      meanratio <- meanratio_paired
   } else {
      meanratio <- meanratio_unpaired
   }
   boot.out <- boot(data, meanratio,R=R, stype="i", strata=data[,"grp"])
   ## Ratio of means (percentile bootstrap) ##
   if(sd(boot.out$t[!is.na(boot.out$t) & is.finite(boot.out$t)])==0){
      conf.int <- c(NA,NA)
   } else {
      conf.int <- boot.ci(boot.out = boot.out, conf = conf.level, type = c("perc"))$perc[4:5]
   }

   estimate <- 1-mean(epg2)/mean(epg1)

   # approximate CI according to Coles for unpaired situation
   talpha<- qt(1-alpha/2,sum(ni)-2)
   lowerCI <- ifelse(estimate==100, NA, 100*(1-mean(epg2)/mean(epg1)*exp(talpha*sqrt(var(epg1)/ni[1]/mean(epg1)^2+var(epg2)/ni[2]/mean(epg2)^2))) )
   upperCI <- ifelse(estimate==100, NA, 100*(1-mean(epg2)/mean(epg1)*exp(-talpha*sqrt(var(epg1)/ni[1]/mean(epg1)^2+var(epg2)/ni[2]/mean(epg2)^2))) )

   return(list(estimate = estimate*100, bootCI = conf.int*100, approxCI=c(lowerCI,upperCI)))
}
