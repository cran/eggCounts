data {
  int J; // number of animals
  array[J] int ystarraw; // McMaster count
  array[J] real CF;
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
  mu ~ gamma(1,0.001);    // priors
  kappa ~ gamma(1,0.7);
  phi ~ beta(1,1);
  mui ~ gamma(kappa, kappa/mu); 
  for (n in 1:J) {             // likelihoods
    if (ystarraw[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | phi), bernoulli_lpmf(0 | phi)+ poisson_lpmf(ystarraw[n] | lambda[n]));
    else
      target += bernoulli_lpmf(0 | phi) + poisson_lpmf(ystarraw[n] | lambda[n]);
  }
}
