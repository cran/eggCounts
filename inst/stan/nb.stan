data {
  int J; // number of animals
  array[J] int ystarraw; // McMaster count
  array[J] real CF;
}
parameters{
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
model{
  mu ~ gamma(1, 0.001);    // priors
  kappa ~ gamma(1, 0.7);
  mui ~ gamma(kappa, kappa/mu);       // likelihoods, gamma(shape, rate)
  ystarraw ~ poisson(lambda);
}

