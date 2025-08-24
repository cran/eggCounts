data{
  int J; // number of animals
  array[J] int ystararaw; // after treatment McMaster count
  array[J] int ystarbraw; // before treatment McMaster count
  array[J] real fpre;
  array[J] real fpost;
}

parameters{
  real<lower=0,upper=1> delta;
  real<lower=0> mu;
}

transformed parameters{
  array[J] real lambdaa;
  array[J] real lambdab;
  for (i in 1:J){
    lambdab[i] = mu/fpre[i];
    lambdaa[i] = delta*mu/fpost[i];
    }
}

model{
  mu ~ gamma(1, 0.001);    // prior
  delta ~ beta(1,1);
  ystararaw ~ poisson(lambdaa);         // likelihoods
  ystarbraw ~ poisson(lambdab);
}
