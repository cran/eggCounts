fecr_probability <- function(stanFit,threshold=0.95,lessthan=TRUE){
  modelName <- stanFit@model_name
  if (modelName %in% c("Bayesian model without zero-inflation","nb","Zero-inflated Bayesian model","zinb")) stop(
    "There is no reduction parameter for this model.")
  fecr<-1-extract(stanFit,"delta")[[1]]
  prob95 <- sum(fecr<threshold)/length(fecr)
  cat(paste0(c("The probability that the reduction is less than",threshold,"is",if(lessthan){prob95*100}else{100-prob95*100},"%.")))
  return(invisible(if(lessthan){prob95*100}else{100-prob95*100}))
}