##' This package implements hierarchical models for faecal egg count data to
##' assess anthelmintic efficacy. Bayesian inference is done via MCMC sampling.
##'
##' \tabular{ll}{
##' Package: \tab eggCounts\cr
##' Type: \tab Package\cr
##' Version: \tab 0.4-1\cr
##' Date: \tab 2015-02-07\cr
##' License: \tab GPL (>= 2) \cr
##' LazyLoad: \tab yes\cr
##' }
##'
# Further details to follow...
##' See \code{demo("fecm", package="eggCounts")} for an illustration.
##' 
##' @name eggCounts-package
##' @aliases eggCounts
##' @docType package
##' @title Hierarchical modelling of faecal egg counts
##' @author Michaela Paul \email{michaela.paul@@uzh.ch}
##' @keywords package
##' @importFrom coda mcmc
##' @importFrom actuar rinvgamma dinvgamma
##' @importFrom boot boot boot.ci
# @importFrom methods setClass setOldClass setGeneric setMethod representation
# signature prototype initialize new
NULL

##' Faecal egg count samples (before and after treatment) 
##'
##' This is an example data set containing 24 eggs per gram (epg) values
##' before and after anthelmintic treatment.
##' The correction factor of the diagnostic technique was 10.
##'
# details...
##'
##' @name epgs
##' @usage data(epgs)
##' @docType data
##' @keywords datasets
NULL

##' Faecal egg count sample 
##'
##' This is an example data set containing 24 eggs per gram (epg) values 
##' of Taenia parasites (Echinococcus) in dogs.
##' The correction factor of the diagnostic technique was 50.
##'
# details...
##'
##' @name echinococcus
##' @usage data(echinococcus)
##' @docType data
##' @keywords datasets
##' @examples
##' data(echinococcus)
##' table(echinococcus)
NULL

##' Abundance of trichostrongyloid eggs in sheep faeces
##' 
##' This data set contains information about the abundance and distribution of 
##' trichostrongyloid eggs in the faeces of 14 groups of commercially farmed 
##' sheep given in Table 1 in Morgan et al. (2005).
##' The faecal egg counts were assumed to follow a negative binomial distribution
##' with mean \eqn{m} and overdispersion parameter \eqn{k}.
##'
##' The data set has columns:
##'
##' \tabular{ll}{
##'   group \tab ID number for the groups \cr
##'   ageclass \tab age class of sheep: "Lambs" or "Ewes" \cr
##'   month \tab month when samples were taken \cr
##'   n \tab number of sheep in group \cr
##'   meanFEC \tab mean number of eggs per gram (epg) of  faeces \cr
##'   k \tab estimated overdispersion parameter \eqn{k} \cr
##'   k.low \tab lower limit of a 95\% confidence interval for \eqn{k} \cr
##'   k.up \tab lower limit of a 95\% confidence interval for \eqn{k} \cr
##'   maxFEC \tab maximal number of eggs per gram of faeces per group \cr
##'   percentageLarger1000 \tab percentage of samples with more than 1000 epg \cr
##'   Chi2 \tab goodness-of-fit statistic for the negative binomial distribution \cr
##'   df \tab corresponding degrees of freedom \cr
##'   p \tab corresponding p-value
##' }
##'
##' @name tab1morgan
##' @usage data(tab1morgan)
##' @format A data frame with 14 rows and 13 variables
##' @source Morgan ER, Cavill L, Curry GE, Wood RM, Mitchell ESE (2005).
##'  Effects of aggregation and sample size on composite faecal egg counts in sheep,
##'  Veterinary Parasitology, 121:79-87.
##' @docType data
##' @keywords datasets
##' @examples
##' data(tab1morgan)
##' if (require("lattice"))
##'    xyplot(k.low+k.up+k ~meanFEC, type="p", pch=19, col=c(8,8,1), data=tab1morgan)
NULL
