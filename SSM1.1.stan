// Focus on the time series of only one firm
// The observations are reported accounting figures, say Gross Profit, (scaled by Rev, total revenue) a vector r of length N (years).
// The covariates are a matrix capital X with K column vectors, incl a unit vector, capturing the  
//   dis/incentives of misreporting the underlying unbiased figures (y) as the observed reported firgures
data {
  int<lower=0> J; // number of firms
  int<lower=0> N; // number of observations (for different years)
  int<lower=0> K; // number of coefficents for the column vectors in covariate matrix X
  int<lower=0> L; // length of simplex vector representing reversal lasting L periods
  real<lower=0> sd_gamma_init;
  real<lower=0> sd_y_init;
  real<lower=0> sd_m_init;
  real alpha_init; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
  real beta_init; // slope coefficient of the AR(1) process of y[n]
  vector[K] g_init;
#  real<lower=0> sd_c_init;
  vector[N] r[J]; // reported figures (scaled by Rev, total revenue)
  matrix[N,K] X[J]; // covariates (capturing the dis/incentives of misreporting)

// forecasts  
//  int<lower=0> N_new; // number of predictions
//  matrix[N_new,K] X_new; // 
}
transformed data {
//  real <lower=0,upper=1> a = 0; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
// GP data
}
parameters {
//  vector[N] y; // underlying unbiased figure (scaled by Rev, total revenue)
  real alpha; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
  real<lower=0,upper=1> beta; // slope coefficient of the AR(1) process of y[n]
  vector[K] g; // coefficients of the K covariates in matrix X
  simplex[L] d; // fractional reversal of prior-period manipluation by accruals
//  real <lower=0,upper=1> d; // fractional reversal of prior-period manipluation by accruals
  real<lower=0> sd_y; // sd of the underlying unbiased figures (vector y)
  real<lower=0> sd_m; // sd of the misreporting extents (vector m = r - y)
  real<lower=0> sd_gamma; // sd of the hyperprior for gamma
  vector[K] err_gamma[J];
//  vector[N] err_y; // underlying unbiased figure (scaled by Rev, total revenue)
  vector[N] err_m[J]; // underlying unbiased figure (scaled by Rev, total revenue)
// real<lower=0,upper=1> integrity;     // integrity level of the society affecting the decision to misreport or not
//  vector[N] integrity;     // integrity level of each CEO affecting the decision to misreport or not
}
transformed parameters {
  vector[N] y[J];
  vector[N] m[J]; // misreporting extent in the reported figures 
  vector[N] D[J];

  for (j in 1:J) {  

// shouldn't good governance restrict |m| or m^2 ? Then impact of X is on sd_m
    m[j] = X[j]*(g + sd_gamma*err_gamma[j]) + sd_m*err_m[j];


    D[j,1] = m[j,1];
    D[j,2] = m[j,2] - d[1]*m[j,1];     
    D[j,3] = m[j,3] - d[1]*m[j,2] - d[2]*m[j,1];  
//    for (l in 1:L) {
//      real reversal = 0;
//      for (i in 1:(l-1))
//        reversal += d[i] * m[j,l-i];
//      D[j,l] = m[j,l] - reversal;
//    }
    
// L = length of simplex vector representing reversal lasting L periods
    D[j,(L+1):N] = m[j,(L+1):N] - d[1]*m[j,L:(N-1)] - d[2]*m[j,(L-1):(N-2)] - d[3]*m[j,(L-2):(N-3)];     // if estimated d < 1, this'd suggest taking more than one-period to reverse
//    for (n in (L+1):N) // L = length of simplex vector representing reversal lasting L periods
//      D[j,n] = m[j,n] - d[1]*m[j,n-1] - d[2]*m[j,n-2] - d[3]*m[j,n-3];     // if estimated d < 1, this'd suggest taking more than one-period to reverse

    y[j] = r[j] - D[j];

  }
//  real ilogodds;     // log odds of the society integrity level
//  ilogodds = logit(integrity);   
}
model {
//  vector[N] m; // misreporting extent in the reported figures 

// priors 
	sd_gamma ~ inv_gamma(2, sd_gamma_init);  // Inverse-gamma (alpha, beta) has mean = beta/(alpha - 1) if alpha > 1
	sd_y ~ inv_gamma(2, sd_y_init);  // Inverse-gamma (alpha, beta) has mean = beta/(alpha - 1) if alpha > 1
	sd_m ~ inv_gamma(2, sd_m_init);  // Inverse-gamma (alpha, beta) has mean = beta/(alpha - 1) if alpha > 1
  alpha ~ normal(alpha_init, 1);
  beta ~ beta(5*beta_init, 5*(1-beta_init));
  g ~ normal(g_init, 1);
  
  for (j in 1:J) { 

    err_gamma[j] ~ normal(0, 1);
    err_m[j] ~ normal(0, 1);   // y should be nonnegative for Rev and nonpositive for Costs

    y[j,2:N] ~ normal(alpha + beta * y[j,1:(N - 1)], sd_y);

  }
//  M ~ bernoulli_logit(ilogodds); 
}
generated quantities {

//  vector[N_new] r_new;
//  for (n in 1:N_new)
//    r_new[n] = normal_rng(Z_new[n] * b, sd_r);
}

