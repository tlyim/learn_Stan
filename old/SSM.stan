// Focus on the time series of only one firm
// The observations are reported accounting figures, say Gross Profit, (scaled by Rev, total revenue) a vector r of length N (years).
// The covariates are a matrix capital X with K column vectors, incl a unit vector, capturing the  
//   dis/incentives of misreporting the underlying unbiased figures (y) as the observed reported firgures
data {
  int<lower=0> N; // number of observations (for different years)
  int<lower=0> K; // number of coefficents for the column vectors in covariate matrix X
  vector[N] r; // reported figures (scaled by Rev, total revenue)
  matrix[N,K] X; // covariates (capturing the dis/incentives of misreporting)

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
  real <lower=0,upper=1> a; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
  real <lower=0,upper=1> b; // slope coefficient of the AR(1) process of y[n]
  vector[K] c; // coefficients of the K covariates in matrix X
  simplex[3] d; // fractional reversal of prior-period manipluation by accruals
//  real <lower=0,upper=1> d; // fractional reversal of prior-period manipluation by accruals
  real<lower=0> sd_y; // sd of the underlying unbiased figures (vector y)
  real<lower=0> sd_m; // sd of the misreporting extents (vector m = r - y)
//  vector[N] err_y; // underlying unbiased figure (scaled by Rev, total revenue)
  vector[N] err_m; // underlying unbiased figure (scaled by Rev, total revenue)
// real<lower=0,upper=1> integrity;     // integrity level of the society affecting the decision to misreport or not
//  vector[N] integrity;     // integrity level of each CEO affecting the decision to misreport or not
}
transformed parameters {
  vector[N] y;
  vector[N] m; // misreporting extent in the reported figures 
  vector[N] D;
  
  m = X*c + err_m;
  D[1] = m[1];
  D[2] = m[2] - d[1]*m[1];     
  D[3] = m[3] - d[1]*m[2] - d[2]*m[1];     
  for (n in 4:N) 
    D[n] = m[n] - d[1]*m[n-1] - d[2]*m[n-2] - d[3]*m[n-3];     // if estimated d < 1, this'd suggest taking more than one-period to reverse
  y[1] = ( a/(1-b) + (r[1] - D[1]) )/2;  
//  y[1] = a/(1-b) + err_y[1];  
  for (n in 2:N) {
    y[n] = ( a + b * y[n-1] + (r[n] - D[n]) )/2;
//    y[n] = a + b * y[n-1] + err_y[n];
//    D[n] = m[n] - d * m[n-1];     // if estimated d < 1, this'd suggest taking more than one-period to reverse
  }


//  m = r - y;

//  real ilogodds;     // log odds of the society integrity level
//  ilogodds = logit(integrity);   
}
model {
//  vector[N] m; // misreporting extent in the reported figures 
  r ~ normal(y + D, sd_y);  // shouldn't good governance restrict |m| or m^2 ? Then impact of X is on sd_m

  err_m ~ normal(0, sd_m);   // y should be nonnegative for Rev and nonpositive for Costs

//  err_y ~ normal(0, sd_y);   // y should be nonnegative for Rev and nonpositive for Costs
//  for (n in 2:N)
//    y[n] ~ normal(alpha + beta * y[n-1], sigma);
//  y[2:N] ~ normal(a + b * y[1:(N - 1)], sd_y);

// c ~ cauchy(0, 2.5); // common prior for each b[K]

// m ~ normal(X*c, sd_m);  // shouldn't good governance restrict |m| or m^2 ? Then impact of X is on sd_m


//  M ~ bernoulli_logit(ilogodds); 

}
generated quantities {

//  vector[N_new] r_new;
//  for (n in 1:N_new)
//    r_new[n] = normal_rng(Z_new[n] * b, sd_r);
}

