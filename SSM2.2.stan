// Focus on the time series of only one firm
// The observations are reported accounting figures, say Gross Profit, (scaled by Rev, total revenue) a vector r of length N (years).
// The covariates are a matrix capital X with K column vectors, incl a unit vector, capturing the  
//   dis/incentives of misreporting the underlying unbiased figures (y) as the observed reported firgures
data {
  int<lower=0> J; // number of firms
  int<lower=0> N; // number of observations (for different years)
  int<lower=0> K; // number of coefficents for the column vectors in covariate matrix X
  int<lower=0> H; // number of coefficents for the column vectors in covariate matrix G
  int<lower=0> L; // length of simplex vector representing reversal lasting L periods
  real<lower=0> sd_gamma_init;
  real<lower=0> sd_omega_init;
  real<lower=0> sd_base_init;
  
//real<lower=0> sd_temp_init;  // for debugging only
  
  real<lower=0> sd_y_init;
//  real<lower=0> sd_m_init;
  real y_init;
  real alpha_init; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
  real beta_init; // slope coefficient of the AR(1) process of y[n]
  vector[K] g_init;
  vector[H] w_init;//  real<lower=0> sd_c_init;

//vector[N] temp[J];  // for debugging only

  vector[N] r[J]; // reported figures (scaled by Rev, total revenue)
  matrix[N,K] X[J]; // covariates (capturing the dis/incentives of misreporting)
  matrix[N,H] G[J]; // covariates (capturing the strength of internal governance and external monitoring mechansims)
//  real<lower=0> base;
// forecasts  
//  int<lower=0> N_new; // number of predictions
//  matrix[N_new,K] X_new; // 
}
transformed data {
real alpha01_init = 0.5*(1 + alpha_init);
//  real <lower=0,upper=1> a = 0; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
}
parameters {
//  vector[N] y; // underlying unbiased figure (scaled by Rev, total revenue)
  real<lower=0,upper=1> alpha01; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
//  real<lower=-1,upper=1> alpha; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
  real<lower=0,upper=1> beta; // slope coefficient of the AR(1) process of y[n]
  vector[K] g; // coefficients of the K covariates in matrix X
  vector[H] w; // coefficients of the H covariates in matrix G
  simplex[L] d; // fractional reversal of prior-period manipluation by accruals
//  real <lower=0,upper=1> d; // fractional reversal of prior-period manipluation by accruals
  real<lower=0> sd_y; // sd of the underlying unbiased figures (vector y)
  real<lower=0> sd_gamma; // sd of the hyperprior for gamma
  real<lower=0> sd_omega; // sd of the hyperprior for omega
real<lower=0> sd_base; // sd of the hyperprior for gamma
real<lower=0,upper=1> mu_base[J]; // sd of the hyperprior for gamma
//vector<lower=0>[N] base[J];  // ,upper=1
vector[N] err_log_base[J];

//real<lower=0> sd_temp; // for debugging only

//  real<lower=0,upper=1> base[J];
  vector[H] err_omega[J];
  vector[K] err_gamma[J];
// vector[K] err_base[J];

//  vector[N] err_m[J]; // underlying unbiased figure (scaled by Rev, total revenue)

// real<lower=0,upper=1> integrity;     // integrity level of the society affecting the decision to misreport or not
//  vector[N] integrity;     // integrity level of each CEO affecting the decision to misreport or not
}
transformed parameters {
  vector[N] y[J];
  vector[N] m[J]; // misreporting extent in the reported figures 
  vector[N] D[J];
//vector<lower=0>[N] base[J];  // ,upper=1
  real<lower=-1,upper=1> alpha; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
  
alpha = 2*alpha01 - 1;
//vector[N] temp0[J];  // for debugging only

  for (j in 1:J) {  
    vector[N] tau[J]; // temptation to misreport
vector[N] rho[J]; // 
vector[N] base[J];  // ,upper=1
//    vector[N] b[J]; // bias effort driven by the temptation to misreport //<lower=-1,upper=1>
//    vector[N] P[J];  // potential room of manipulation constrained by governance mechanisms //<lower=0>

// shouldn't good governance restrict |m| or m^2 ? Then impact of X is on sd_m
//    m[j] = X[j]*(g + sd_gamma*err_gamma[j]);  // + sd_m*err_m[j];
    tau[j] = X[j]*(g + sd_gamma*err_gamma[j]);  // temptation to misreport

    rho[j] = G[j]*(w + sd_omega*err_omega[j]);

base[j] = exp( mu_base[j] + sd_base*err_log_base[j] );

 //    b[j] = (exp(tau[j]) - 1) ./ (exp(tau[j]) + 1);
//    P[j] = inv_logit( (-1)*G[j]*(w + sd_omega*err_omega[j]) );

//    m[j] = b[j] .* P[j] .* rep_vector(base[j], N);  // extent of misreporting resulting from bias effort 
    m[j] = ( base[j] .* (exp(tau[j]) - 1) ) ./ ( 1 + exp(tau[j]) + exp(rho[j]) + exp(tau[j]+rho[j]) );

    
//    D[j,1] = m[j,1];
//    D[j,2] = m[j,2] - d[1]*m[j,1];     
//    D[j,3] = m[j,3] - d[1]*m[j,2] - d[2]*m[j,1];  
    for (l in 1:L) {
      real reversal = 0;
      for (i in 1:(l-1))
        reversal += d[i] * m[j,l-i];
      D[j,l] = m[j,l] - reversal;
      }
// L = length of simplex vector representing reversal lasting L periods
//    D[j,(L+1):N] = m[j,(L+1):N] - d[1]*m[j,L:(N-1)] - d[2]*m[j,(L-1):(N-2)] - d[3]*m[j,(L-2):(N-3)];     // if estimated d < 1, this'd suggest taking more than one-period to reverse
    D[j,(L+1):N] = m[j,(L+1):N];
    for (l in 1:L) 
//      D[j,(L+1):N] += (-1)* d[l] * m[j,(L+1-l):(N-l)];
      D[j,(L+1):N] = D[j,(L+1):N] + (-1)* d[l] * m[j,(L+1-l):(N-l)];


    y[j] = r[j] - D[j];

  }

//  real ilogodds;     // log odds of the society integrity level
//  ilogodds = logit(integrity);   
}
model {

//vector[N] log_base[J];

// priors 
	sd_y ~ normal(sd_y_init, 0.12);  // Inverse-gamma (alpha, beta) has mean = beta/(alpha - 1) if alpha > 1
//	sd_y ~ inv_gamma(2, sd_y_init);  // Inverse-gamma (alpha, beta) has mean = beta/(alpha - 1) if alpha > 1
	sd_gamma ~ normal(sd_gamma_init, 0.025);  // Inverse-gamma (alpha, beta) has mean = beta/(alpha - 1) if alpha > 1
//	sd_omega ~ inv_gamma(2, sd_omega_init);  // Inverse-gamma (alpha, beta) has mean = beta/(alpha - 1) if alpha > 1
	sd_omega ~ normal(sd_omega_init, 0.045);  // Inverse-gamma (alpha, beta) has mean = beta/(alpha - 1) if alpha > 1
sd_base ~ normal(sd_base_init, 0.385); //0.2);  // Inverse-gamma (alpha, beta) has mean = beta/(alpha - 1) if alpha > 1

//  base ~ inv_gamma(2, y_init);    // for debugging only
//sd_temp ~ inv_gamma(2, sd_temp_init);    // for debugging only

//  alpha ~ normal(alpha_init, 0.1);
  alpha01 ~ beta(3*alpha01_init, 3*(1-alpha01_init));
  beta ~ beta(6*beta_init, 6*(1-beta_init));
  g ~ normal(g_init, 0.2);
  w ~ normal(w_init, 0.2); //1);
  
//  d ~ dirichlet(rep_vector(0.1, L));   // =1.0 means the dist is uniform over the possible simplices;<1.0 toward corners 
  //  V ~ uniform(-1, 1);

  
  for (j in 1:J) { 

    err_gamma[j] ~ normal(0, 1);
    err_omega[j] ~ normal(0, 1);
//    err_base[j] ~ normal(0, 1);

  err_log_base[j] ~ normal(0, 1);

//    log(base[j]) ~ normal(mu_base[j], sd_base);
//    target += -log(base[j]);
    
    y[j,2:N] ~ normal(alpha + beta * y[j,1:(N - 1)], sd_y);

  }
//  M ~ bernoulli_logit(ilogodds); 
}
generated quantities {

//  vector[N_new] r_new;
//  for (n in 1:N_new)
//    r_new[n] = normal_rng(Z_new[n] * b, sd_r);
}

