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
//  vector[K] gamma[J];  //gamma is the 'random cofficient c' normally distributed with mean vector c and sd vector sd_gamma
//  vector[K] err_gamma[J];


  for (j in 1:J) {
    vector[K] err_gamma[J];
    vector[N] P[J];  //<lower=0> potential room of manipulation constrained by governance mechanisms
    vector[N] tau[J]; // temptation to misreport
    vector[N] b[J]; //<lower=-1,upper=1> bias effort driven by the temptation to misreport
    
    for (k in 1:K) {
      err_gamma[j,k] = normal_rng(0, 1);
    }
    
    for (n in 1:N){
      P[j,n] = exp(2); 
// shouldn't good governance restrict |m| or m^2 ? Then impact of X is on sd_m
//    m[j] = X[j]*(g + sd_gamma*err_gamma[j]);  // + sd_m*err_m[j];
      tau[j,n] = X[j,n]*(g + sd_gamma*err_gamma[j]);  // temptation to misreport
    }

    b[j] = (exp(tau[j]) - 1) ./ (exp(tau[j]) + 1);

//    m[j,1] = normal_rng(X[j,1]*gamma[j], sd_m);
//    m[j,1] = X[j,1]*gamma[j];
    m[j,1] = b[j,1] .* P[j,1];  // extent of misreporting resulting from bias effort and 
    
    y[j,1] = normal_rng(alpha/(1-beta), sd_y);   // y should be nonnegative for Rev and nonpositive for Costs
// Perhaps best to model y as always positive, with IsCost = 1, 0 to indicate Cost or Rev item 
    for (n in 2:N) {
//      m[j,n] = normal_rng(X[j,n]*gamma[j], sd_m);
//      m[j,n] = X[j,n]*gamma[j];
      m[j,n] = b[j,n] .* P[j,n];  // extent of misreporting resulting from bias effort and 

      y[j,n] = alpha + normal_rng(beta*y[j,n-1], sd_y);
    }
  }

}


