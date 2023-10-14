data{
  int Ja; // number of animals
  int Jb;
  array[Ja] int ystararaw; // after treatment McMaster count
  array[Jb] int ystarbraw; // before treatment McMaster count
  array[Ja] real fpost;
  array[Jb] real fpre;
}
parameters{
  real<lower=0> kappa;
  real<lower=0> mu;
  real<lower=0,upper=1> delta;
  array[Jb] real<lower=0> mub; // true epg before treatment
  array[Ja] real<lower=0> mua; // true epg after treatment
}
transformed parameters{
  array[Ja] real lambdaa;
  array[Jb] real lambdab;
  for (i in 1:Jb){
    lambdab[i] = mub[i]/fpre[i];
  }
  for (i in 1:Ja){
    lambdaa[i] = delta*mua[i]/fpost[i];
  }
}
model{
  mu ~ gamma(1,0.001);    // priors
  kappa ~ gamma(1,0.7);
  delta ~ beta(1,1);
  mub ~ gamma(kappa, kappa/mu); // likelihoods
  mua ~ gamma(kappa, kappa/mu);   
  ystarbraw ~ poisson(lambdab);
  ystararaw ~ poisson(lambdaa);
}
