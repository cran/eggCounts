# set default values for the priors 
fec_setPrior <- function(muPrior, kappaPrior, phiPrior){
  if(missing(muPrior)) muPrior = list(priorDist = "gamma",hyperpars=c(1,0.001))
  if(is.null(muPrior[["priorDist", exact = TRUE]])) muPrior$priorDist = "gamma"
  if(is.null(muPrior[["hyperpars", exact = TRUE]])) muPrior$hyperpars = c(1,0.001)
  
  if(missing(kappaPrior)) kappaPrior = list(priorDist = "gamma",hyperpars=c(1,0.7))
  if(is.null(kappaPrior[["priorDist", exact = TRUE]])) kappaPrior$priorDist = "gamma"
  if(is.null(kappaPrior[["hyperpars", exact = TRUE]])) kappaPrior$hyperpars = c(1,0.7)
  
  if(missing(phiPrior)) phiPrior = list(priorDist = "beta",hyperpars=c(1,1))
  if(is.null(phiPrior[["priorDist", exact = TRUE]])) phiPrior$priorDist = "beta"
  if(is.null(phiPrior[["hyperpars", exact = TRUE]])) phiPrior$hyperpars = c(1,1)
  
  return(list(mu=muPrior,kappa=kappaPrior,phi = phiPrior))
}

# Stan model code for without zero inflation  -----------------------------

nb_stan <- function(priors){
  #hyperparameters for pre-treatment mean mu
  a.mu <- priors$mu$hyperpars[1]
  b.mu <- priors$mu$hyperpars[2]
  dist.mu <- priors$mu$priorDist
  #hyperparameters  for overdispersion parameter kappa
  a.kappa <- priors$kappa$hyperpars[1]
  b.kappa <- priors$kappa$hyperpars[2]
  dist.kappa <- priors$kappa$priorDist
  paste0('data {
         int J; // number of animals
         array[J] int ystarraw; // McMaster count
         array[J] int CF;
         }
         parameters {
          real<lower=0> kappa;
          real<lower=0> mu;
          array[J] real<lower=0> mui;
         }
         transformed parameters{
          array[J] real lambda;
          for (i in 1:J){
          lambda[i] = mui[i]/CF[i];
          }
         }
         model {
          mu ~ ',dist.mu,'(',a.mu,',',b.mu,');    // prior
          kappa ~ ',dist.kappa,'(',a.kappa,',',b.kappa,');
          mui ~ gamma(kappa,kappa/mu);       // likelihoods
          ystarraw ~ poisson(lambda);
        }')
}

# Stan model code for with zero inflation  --------------------------------

zinb_stan <- function(priors){
  #hyperparameters for pre-treatment mean mu
  a.mu <- priors$mu$hyperpars[1]
  b.mu <- priors$mu$hyperpars[2]
  dist.mu <- priors$mu$priorDist
  #hyperparameters  for overdispersion parameter kappa
  a.kappa <- priors$kappa$hyperpars[1]
  b.kappa <- priors$kappa$hyperpars[2]
  dist.kappa <- priors$kappa$priorDist
  #hyperparameters for prevalence 
  a.phi <- priors$phi$hyperpars[1]
  b.phi <- priors$phi$hyperpars[2]
  dist.phi <- priors$phi$priorDist
  paste0('data {
         int J; // number of animals
         array[J] int ystarraw; // McMaster count
         array[J] int CF;
         }
         parameters {
          real<lower=0> kappa;
          real<lower=0> mu;
          array[J] real<lower=0> mui;
          real<lower=0,upper=1> phi;
         }
         transformed parameters{
          array[J] real lambda;
          for (i in 1:J){
          lambda[i] = mui[i]/CF[i];
          }
         }
          model {
          mu ~ ',dist.mu,'(',a.mu,',',b.mu,');       // prior
          kappa ~ ',dist.kappa,'(',a.kappa,',',b.kappa,');
          phi ~ ',dist.phi,'(',a.phi,',',b.phi,');
          mui ~ gamma(kappa, kappa/mu); 
          for (n in 1:J) {             // likelihoods
          if (ystarraw[n] == 0)
          target += log_sum_exp(bernoulli_lpmf(1 | phi), bernoulli_lpmf(0 | phi)+poisson_lpmf(ystarraw[n] | lambda[n]));
          else
          target += bernoulli_lpmf(0 | phi) + poisson_lpmf(ystarraw[n] | lambda[n]);
          }
          }')
}