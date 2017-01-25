data {
  int J; // number of animals
  int ystararaw[J]; // after treatment McMaster count
  int ystarbraw[J]; // before treatment McMaster count
  int fpre[J];
  int fpost[J];
}
parameters {
  real<lower=0> kappa;
  real<lower=0> mu;
  real<lower=0,upper=1> delta;
  real<lower=0> mub[J];
}
transformed parameters{
  real lambdaa[J];
  real lambdab[J];
  real kappamu;
  for (i in 1:J){
    lambdab[i] = mub[i]/fpre[i];
  }
  for (i in 1:J){
    lambdaa[i] = delta*mub[i]/fpost[i];
  }  
  kappamu = kappa/mu;
}
model {
  mu ~ gamma(1,0.001);    // prior
  kappa ~ gamma(1,0.7);
  delta ~ beta(1,1);
  target += gamma_lpdf(mub | kappa,kappamu);   // likelihoods 
  ystararaw ~ poisson(lambdaa);
  ystarbraw ~ poisson(lambdab);

}
