data {
  int J; // number of animals
  int ystarraw[J]; // McMaster count
  int CF[J];
}
parameters {
  real<lower=0> kappa;
  real<lower=0> mu;
  real<lower=0> mui[J];
  real<lower=0,upper=1> phi;
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
  phi ~ beta(1,1);
  mui ~ gamma(kappa,kappamu); 
  for (n in 1:J) {             // likelihoods
    if (ystarraw[n] == 0)
      increment_log_prob(log_sum_exp(bernoulli_log(1,phi), bernoulli_log(0,phi)+poisson_log(ystarraw[n],lambda[n])));
    else
      increment_log_prob(bernoulli_log(0,phi) + poisson_log(ystarraw[n],lambda[n]));
  }
}
