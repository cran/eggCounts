stan2mcmc <- function(stanFit){
  modelName <- stanFit@model_name
  switch(modelName,
         "Zero-inflated Bayesian model for paired design"=,
         "zipaired"={
           meanEPG.untreated<-extract(stanFit,"mu")[[1]]*(1-extract(stanFit,"phi")$phi)
           meanEPG.treated<-extract(stanFit,"mu")[[1]]*extract(stanFit,"delta")$delta*(1-extract(stanFit,"phi")$phi)
           fecr<-1-extract(stanFit,"delta")[[1]]
           result<-cbind(fecr,meanEPG.untreated,meanEPG.treated)
           output<-cbind(result,as.data.frame(extract(stanFit,c("kappa","phi","delta"))))
         },
         "Zero-inflated Bayesian model for unpaired design"=,
         "ziunpaired"={
           meanEPG.untreated<-extract(stanFit,"mu")[[1]]*(1-extract(stanFit,"phi")$phi)
           meanEPG.treated<-extract(stanFit,"mu")[[1]]*extract(stanFit,"delta")$delta*(1-extract(stanFit,"phi")$phi)
           fecr<-1-extract(stanFit,"delta")[[1]]
           result<-cbind(fecr,meanEPG.untreated,meanEPG.treated)
           output<-cbind(result,as.data.frame(extract(stanFit,c("kappa","phi","delta"))))
         },
         "Bayesian model without zero-inflation for paired design"=,
         "paired"={
           meanEPG.untreated<-extract(stanFit,"mu")[[1]]
           meanEPG.treated<-extract(stanFit,"mu")[[1]]*extract(stanFit,"delta")$delta
           fecr<-1-extract(stanFit,"delta")[[1]]
           result<-cbind(fecr,meanEPG.untreated,meanEPG.treated)
           output<-cbind(result,as.data.frame(extract(stanFit,c("kappa","delta"))))
         },
         "Bayesian model without zero-inflation for unpaired design"=,
         "unpaired"={
           meanEPG.untreated<-extract(stanFit,"mu")[[1]]
           meanEPG.treated<-extract(stanFit,"mu")[[1]]*extract(stanFit,"delta")$delta
           fecr<-1-extract(stanFit,"delta")[[1]]
           result<-cbind(fecr,meanEPG.untreated,meanEPG.treated)
           output<-cbind(result,as.data.frame(extract(stanFit,c("kappa","delta"))))
         },
         "Bayesian model without zero-inflation"=,
         "nb"={
           meanEPG<-extract(stanFit,"mu")[[1]]
           kappa<-extract(stanFit,"kappa")$kappa
           output<-cbind(meanEPG=meanEPG,kappa=kappa)
         },
         "Zero-inflated Bayesian model"=,
         "zinb"={
           meanEPG<-extract(stanFit,"mu")[[1]]*(1-extract(stanFit,"phi")[[1]])
           phi<-extract(stanFit,"phi")$phi
           kappa<-extract(stanFit,"kappa")$kappa
           output<-cbind(meanEPG=meanEPG,kappa=kappa,phi=phi)
         }
        )
  return(invisible(mcmc(output,thin=stanFit@sim$thin)))
}