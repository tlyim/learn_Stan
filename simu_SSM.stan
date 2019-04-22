data {
  int<lower=0> N; // number of observations (for different years)
  int<lower=0> K; // number of coefficents for the column vectors in covariate matrix X
  matrix[N,K] X; // covariates (capturing the dis/incentives of misreporting)

  real a; // intercept coefficient (drift) of the AR(1) process of the underlying unbiased accounting figure y[n]
  real b; // slope coefficient of the AR(1) process of y[n]
  vector[K] c; // coefficients of the K covariates in matrix X
  real<lower=0> sd_y; // sd of the underlying unbiased figure (vector y)
  real<lower=0> sd_m; // sd of the misreporting extent (vector m = r - y)
}
transformed data {
//  real a = 0;
}
parameters {}
model {}
generated quantities {
  vector[N] y; // underlying unbiased figures (eg, Gross Profit, scaled by Rev)
  vector[N] m; // misreporting extent in the reported figures 

  m[1] = normal_rng(X[1]*c, sd_m); 
  y[1] = normal_rng(a/(1-b), sd_y);   // y should be nonnegative for Rev and nonpositive for Costs
// Perhaps best to model y as always positive, with IsCost = 1, 0 to indicate Cost or Rev item 
  for (n in 2:N) {
    m[n] = normal_rng(X[n]*c, sd_m); 
    y[n] = a + normal_rng(b*y[n-1], sd_y);
  }
}


