data {
  int<lower=0> J; // number of firms
  int<lower=0> N; // number of observations (for different years)
  int<lower=0> K; // number of coefficents for the column vectors in covariate matrix X
  matrix[N,K] X[J]; // covariates (capturing the dis/incentives of misreporting)

  real alpha; // intercept coefficient (drift) of the AR(1) process of the underlying unbiased accounting figure y[n]
  real beta; // slope coefficient of the AR(1) process of y[n]
  vector[K] g; // coefficients of the K covariates in matrix X
  real<lower=0> sd_gamma; // sd of the hyperprior for gamma
  real<lower=0> sd_y; // sd of the underlying unbiased figure (vector y)
//  real<lower=0> sd_m; // sd of the misreporting extent (vector m = r - y)
}
transformed data {
//  real a = 0;
}
parameters {}
model {}
generated quantities {
  vector[N] y[J]; // underlying unbiased figures (eg, Gross Profit, scaled by Rev)
  vector[N] m[J]; // misreporting extent in the reported figures 
  vector[K] gamma[J];  //gamma is the 'random cofficient c' normally distributed with mean vector c and sd vector sd_gamma

  for (j in 1:J) {
    for (k in 1:K) {
      gamma[j,k] = normal_rng(g[k], sd_gamma);
    }
//    m[j,1] = normal_rng(X[j,1]*gamma[j], sd_m);
    m[j,1] = X[j,1]*gamma[j];
    y[j,1] = normal_rng(alpha/(1-beta), sd_y);   // y should be nonnegative for Rev and nonpositive for Costs
// Perhaps best to model y as always positive, with IsCost = 1, 0 to indicate Cost or Rev item 
    for (n in 2:N) {
//      m[j,n] = normal_rng(X[j,n]*gamma[j], sd_m);
      m[j,n] = X[j,n]*gamma[j];
      y[j,n] = alpha + normal_rng(beta*y[j,n-1], sd_y);
    }
  }

}


