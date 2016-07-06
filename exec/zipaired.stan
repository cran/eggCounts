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
  real<lower=0,upper=1> phi;
}
transformed parameters{
  real lambdaa[J];
  real lambdab[J];
  real kappamu;
  for (i in 1:J){
    lambdab[i] = mub[i]/fpre[i];
    lambdaa[i] = delta*mub[i]/fpost[i];
  }
  kappamu = kappa/mu;
}
model {
  mu ~ gamma(1,0.001);    // prior
  kappa ~ gamma(1,0.7);
  delta ~ beta(1,1);
  phi ~ beta(1,1);
  mub ~ gamma(kappa,kappamu);   // likelihoods
  for (n in 1:J) {
    if (ystarbraw[n] == 0)
      target +=  log_sum_exp(bernoulli_lpmf(1 | phi), bernoulli_lpmf(0 | phi)+poisson_lpmf(ystarbraw[n] | lambdab[n]));
    else
      target += bernoulli_lpmf(0 | phi) + poisson_lpmf(ystarbraw[n] | lambdab[n]);
  }
  for (n in 1:J) {
    if (ystararaw[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | phi), bernoulli_lpmf(0 | phi)+poisson_lpmf(ystararaw[n] | lambdaa[n]));
    else 
      target += bernoulli_lpmf(0 | phi) + poisson_lpmf(ystararaw[n] | lambdaa[n]);
  }
}
