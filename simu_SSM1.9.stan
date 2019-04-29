data {
  int<lower=0> J; // number of firms
  int<lower=0> N; // number of observations (for different years)
  int<lower=0> K; // number of coefficents for the column vectors in covariate matrix X
  int<lower=0> H; // number of coefficents for the column vectors in covariate matrix G
  matrix[N,K] X[J]; // covariates (capturing the dis/incentives of misreporting)
  matrix[N,H] G[J]; // covariates (capturing the dis/incentives of misreporting)

  real alpha; // intercept coefficient (drift) of the AR(1) process of the underlying unbiased accounting figure y[n]
  real beta; // slope coefficient of the AR(1) process of y[n]
  vector[K] g; // coefficients of the K covariates in matrix X
  vector[H] w; // coefficients of the H covariates in matrix G

real<lower=0> sd_temp; // for debugging only

  real<lower=0> sd_omega; // sd of the hyperprior for gamma
  real<lower=0> sd_gamma; // sd of the hyperprior for gamma
  real<lower=0> sd_y; // sd of the underlying unbiased figure (vector y)
  real<lower=0> base;
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
  vector[N] b[J]; //<lower=-1,upper=1> bias effort driven by the temptation to misreport
  vector[N] P[J];  //<lower=0> potential room of manipulation constrained by governance mechanisms
vector[N] temp[J];  // for debugging only

  for (j in 1:J) {
    vector[K] err_gamma[J];
    vector[H] err_omega[J];
    vector[N] tau[J]; // temptation to misreport
    
    for (k in 1:K) {
      err_gamma[j,k] = normal_rng(0, 1);
      }
    for (h in 1:H) {
      err_omega[j,h] = normal_rng(0, 1);
      }

    tau[j] = X[j]*(g + sd_gamma*err_gamma[j]);  // temptation to misreport
    b[j] = (exp(tau[j]) - 1) ./ (exp(tau[j]) + 1);
//    P[j] = exp( (-1)*G[j]*(w + sd_omega*err_omega[j]) );
//    P[j] = 1 ./ square( G[j]*(w + sd_omega*err_omega[j]) );
    P[j] = rep_vector(base, N) + inv_square( G[j]*(w + sd_omega*err_omega[j]) );
    //*inv_square( G[j]*(w + sd_omega*err_omega[j]) ) ./ inv_square( G[j]*(w + sd_omega*err_omega[j]) );
temp[j] = normal_rng(0, sd_temp) + G[j]*(w);// + sd_omega*err_omega[j]) ; // for debugging only

//    m[j,1] = normal_rng(X[j,1]*gamma[j], sd_m);
//    m[j,1] = X[j,1]*gamma[j];
    m[j] = b[j] .* P[j];  // extent of misreporting resulting from bias effort and 
//    m[j,1] = b[j,1] .* P[j,1];  // extent of misreporting resulting from bias effort and 
    
    y[j,1] = normal_rng(alpha/(1-beta), sd_y);   // y should be nonnegative for Rev and nonpositive for Costs
// Perhaps best to model y as always positive, with IsCost = 1, 0 to indicate Cost or Rev item 
    for (n in 2:N) {
//      m[j,n] = b[j,n] .* P[j,n];  // extent of misreporting resulting from bias effort and 
      y[j,n] = alpha + normal_rng(beta*y[j,n-1], sd_y);
    }
  }

}


