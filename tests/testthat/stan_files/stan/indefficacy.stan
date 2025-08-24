data{
  int J; // number of animals
  array[J] int ystararaw; // after treatment McMaster count
  array[J] int ystarbraw; // before treatment McMaster count
  array[J] real fpre;
  array[J] real fpost;
}
parameters{
  real<lower=0> kappa;
  real<lower=0> mu;
  array[J] real<lower=0> delta;
  real<lower=0> delta_shape;
  real<lower=0, upper=1> delta_mu;
  array[J] real<lower=0> mub;
}
transformed parameters{
  array[J] real lambdaa;
  array[J] real lambdab;
  for (i in 1:J){
    lambdab[i] = mub[i]/fpre[i];
    lambdaa[i] = delta[i]*mub[i]/fpost[i];
  }
}
model{
  mu ~ gamma(1,0.001);    // prior
  kappa ~ gamma(1,0.7);
  delta ~ gamma(delta_shape, delta_shape/delta_mu); // shape, rate
  delta_shape ~ normal(2, 1);
  delta_mu ~ beta(1,1);
  mub ~ gamma(kappa, kappa/mu);   // likelihoods 
  ystararaw ~ poisson(lambdaa);
  ystarbraw ~ poisson(lambdab);
}
