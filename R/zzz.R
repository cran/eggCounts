 .onLoad <- function(libname, pkgname) { # nocov start
   if (!("methods" %in% .packages())) attachNamespace("methods")
     modules <- paste0("stan_fit4", names(stanmodels), "_mod")
     for (m in modules) loadModule(m, what = TRUE)
   } 
.onAttach <- function(...) {
  rstanarmLib <- dirname(system.file(package = "eggCounts"))
#  pkgdesc <- utils::packageDescription("eggCounts", lib.loc = eggCountsLib)
#  builddate <- gsub(';.*$', '', pkgdesc$Packaged)
#  packageStartupMessage(paste("eggCounts (Version ", pkgdesc$Version, ")", sep = ""))
  packageStartupMessage("- The compiling time for the first time using non-default priors can be up to 20s.")
}

