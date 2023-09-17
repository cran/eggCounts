.onAttach <- function(...){
  packageStartupMessage(paste0("Loading eggCounts (version ", utils::packageVersion("eggCounts"),"):
- The compilation time for the first time using non-default priors can be up to 20s.
- For a graphical user interface of the package implemented in R Shiny, visit: http://shiny.math.uzh.ch/user/furrer/shinyas/shiny-eggCounts/" ))
}

