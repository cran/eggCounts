data {
  int J; // number of animals
  int ystarraw[J]; // McMaster count
  int CF[J];
}
parameters {
  real<lower=0> kappa;
  real<lower=0> mu;
  real<lower=0> mui[J];
}
transformed parameters{
  real lambda[J];
  real kappamu;
  for (i in 1:J){
    lambda[i] <- mui[i]/CF[i];
  }
  kappamu <- kappa/mu;
}
model {
  mu ~ gamma(1,0.001);    // prior
  kappa ~ gamma(1,0.7);
  mui ~ gamma(kappa,kappamu);       // likelihoods
  ystarraw ~ poisson(lambda);
}

