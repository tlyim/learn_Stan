
// RStan provides the expose_stan_functions function for exporting such functions to the R global environment so that they can be tested in R to ensure they are working properly. 

data {
  int<lower=0> J;          // number of schools 
  real y[J];               // estimated treatment effects
  real<lower=0> sigma[J];  // s.e. of effect estimates 
}

parameters {
  real mu; 
  real<lower=0> tau;
  vector[J] eta;
}

transformed parameters {
  vector[J] theta;
  theta = mu + tau * eta;
}

model {
  target += normal_lpdf(eta | 0, 1);
  target += normal_lpdf(y | theta, sigma);
}
