#####################################################
## compute standard FECR, CI based on Levecke 2015
#####################################################

fecrtCI_Levecke <- function(epg1, epg2, paired=FALSE, alpha=0.05, R=1999){
  if (!is.logical(paired))
    stop("The paired argument must be a logical", call.=FALSE)
  if ((R < 1)|(ceiling(R)!=floor(R)))
    stop("'R' should be a positive integer", call.=FALSE)
  if ((alpha<0)||(alpha>1))
    stop("alpha must be between 0 and 1", call.=FALSE)
  if (R > 9999)      cat("NOTE: 'R' seems high\n")
  if (sum(epg2, na.rm = TRUE) == 0) message(cat("NOTE: Confidence interval cannot be computed when all of post-treatment counts are 0. \nPlease use Bayesian models when interval estimates are required.\n"))
  if (any(is.na(c(epg1,epg2)))) message(cat("NA values were removed upon evaluation.\n"))
  
  conf.level <- 1-alpha
  lowerCI <- NA
  upperCI <- NA
  
  if(paired){ # paired design
    nas <- is.na(epg1) | is.na(epg2)
    epg1 <- epg1[!nas]
    epg2 <- epg2[!nas]
    
    ni <- c(length(epg1), length(epg2))
    data <- data.frame(epgs = c(epg1,epg2), grp = as.factor(rep(1:2, ni)))
    
    meanratio_paired <- function(d,i){
      ind1 <- i[1:ni[1]]
      ind2 <- ind1+ni[1]
      x1 <- d[ind1, 1]
      x2 <- d[ind2, 1]
      1-mean(x2)/mean(x1)
    }
    boot.out <- boot(data, meanratio_paired, R=R, stype="i", strata=data[,"grp"])
    ## Ratio of means (percentile bootstrap) ##
    if (sum(epg1)!=0){
      if(sd(boot.out$t[!is.na(boot.out$t) & is.finite(boot.out$t)])==0){
        conf.int <- c(NA,NA)
      } else {
        conf.int <- boot.ci(boot.out = boot.out, conf = conf.level, type = c("perc"))$perc[4:5]
      }
      } else {conf.int <- c(NA,NA)}
    
    estimate <- 1-mean(epg2)/mean(epg1)
    talpha<- qt(1-alpha/2,sum(ni)-2)
    
    if(estimate != 1) {
      estimate.variance <- (mean(epg2)/mean(epg1))^2 * (
        var(epg2)/(ni[2] * mean(epg2)^2) + var(epg1)/(ni[1] * mean(epg1)^2) - 
            2 * mean(sum((epg2 - mean(epg2))*(epg1 - mean(epg1))))/(ni[1] * ni[2] * mean(epg2) * mean(epg1))
        )
      
      lowerCI <- 100 * (estimate - talpha * sqrt(estimate.variance))
      upperCI <- 100 * (estimate + talpha * sqrt(estimate.variance))
    }
  } else { # unpaired design 
    epg1 <- epg1[!is.na(epg1)]
    epg2 <- epg2[!is.na(epg2)]
    
    ni <- c(length(epg1), length(epg2))
    data <- data.frame(epgs = c(epg1,epg2), grp = as.factor(rep(1:2, ni)))
    meanratio_unpaired <- function(d, i) {
      ind1 <- i[1:ni[1]]
      ind2 <- i[-(1:ni[1])]
      x1 <- d[ind1, 1]
      x2 <- d[ind2, 1]
      1-mean(x2)/mean(x1)
    }
    
    boot.out <- boot(data, meanratio_unpaired, R=R, stype="i", strata=data[,"grp"])
    ## Ratio of means (percentile bootstrap) ##
    if (sum(epg1)!=0){
      if(sd(boot.out$t[!is.na(boot.out$t) & is.finite(boot.out$t)])==0){
        conf.int <- c(NA,NA)
      } else {
        conf.int <- boot.ci(boot.out = boot.out, conf = conf.level, type = c("perc"))$perc[4:5]
      }} else {conf.int <- c(NA,NA)}
    
    estimate <- 1-mean(epg2)/mean(epg1)
    talpha<- qt(1-alpha/2,sum(ni)-2)
    
    if(estimate != 1) {
      estimate.variance <- (mean(epg2)/mean(epg1))^2 * (
          var(epg1)/(mean(epg1)^2 * ni[1]) + var(epg2)/(mean(epg2)^2 * ni[2])
          )
      
      lowerCI <- 100 * (estimate - talpha * sqrt(estimate.variance))
      upperCI <- 100 * (estimate + talpha * sqrt(estimate.variance))
    }
  }
  return(list(estimate = estimate*100, bootCI = conf.int*100, approxCI=c(lowerCI,upperCI)))
}



