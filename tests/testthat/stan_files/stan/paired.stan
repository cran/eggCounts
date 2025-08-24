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
  real<lower=0,upper=1> delta;
  array[J] real<lower=0> mub;
}
transformed parameters{
  array[J] real lambdaa;
  array[J] real lambdab;
  for (i in 1:J){
    lambdab[i] = mub[i]/fpre[i];
    lambdaa[i] = delta*mub[i]/fpost[i];
  }
}
model{
  mu ~ gamma(1,0.001);    // priors
  kappa ~ gamma(1,0.7);
  delta ~ beta(1,1);
  mub ~ gamma(kappa, kappa/mu);   // likelihoods, gamma(shape, rate)
  ystararaw ~ poisson(lambdaa);
  ystarbraw ~ poisson(lambdab);

}
