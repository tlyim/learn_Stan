data {
  int<lower=0> J; // number of firms
  int<lower=0> N; // number of observations (for different years)
  int<lower=0> Q; // number of seasons, eg, 4 for quarters
  int<lower=1,upper=Q> q[N]; // index for seasons, eg, q[n] = 1 if n belongs to Q1
  int<lower=0> I; // number of coefficents for the column vectors in covariate matrix Z (which can equal X)
  int<lower=0> K; // number of coefficents for the column vectors in covariate matrix X
  int<lower=0> H; // number of coefficents for the column vectors in covariate matrix G
  int<lower=0> L; // length of simplex vector representing reversal lasting L periods

real rho; // coeff. of the LT impact of real EM
  vector[N] r[J]; // reported figures (scaled by Rev, total revenue)
  matrix[N,I] Z[J]; // covariates (capturing the dis/incentives of real earnings management)
  matrix[N,K] X[J]; // covariates (capturing the dis/incentives of misreporting)
  matrix[N,H] G[J]; // covariates (capturing the strength of internal governance and external monitoring mechansims)
}
transformed data {
}
parameters {
  real mu_u1;
//  real mu_alpha; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
real alpha_ratio_raw; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
  real<lower=0,upper=1> beta; // slope coefficient of the AR(1) process of y[n]
  real<lower=0> sd_y; // sd of the underlying unbiased figures (vector y)
  vector[Q-1] season_raw[J];
  vector[Q-1] mu_season;   // mu of the hyperprior for season_raw
  real<lower=0> sd_season; // sd of the hyperprior for season_raw

vector[I] p; // coefficients of the K covariates in matrix Z
real<lower=0> sd_pi; // sd of the hyperprior for pi
vector[I] err_pi[J];  // error term in p
//real<lower=0,upper=1> rho_raw; // coeff. of the LT impact of real EM
}
transformed parameters {
//vector[N] u[J];
//vector[N] err[J];  // error term in p
  vector[Q] season_q[J];
  vector[N] season_n[J];
vector[N] realLT[J]; // LT impact of real EM on alpha[j,n]
real mu_alpha; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
real alpha_ratio; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]

alpha_ratio = exp(alpha_ratio_raw);
//rho = exp(2*rho_raw);
mu_alpha = rho * alpha_ratio;




  for (j in 1:J) {
    vector[N] zeta[J]; // temptation to manage current-period real earnings upward
//===============================================================================
// Define real EM component
    zeta[j] = Z[j]*(p + sd_pi*err_pi[j]);  // temptation to manage current-period real earnings upward
//    RealLT[j] = log( log1p_exp((-rho)*zeta[j]) - log1p_exp((-rho)*(zeta[j]+1)) ) - log(rho);
//    RealLT[j] = ( log1p_exp((-rho)*zeta[j]) - log1p_exp((-rho)*(zeta[j]+1)) )/rho;
    realLT[j] = ( log1p_exp((-rho)*zeta[j]) - log1p_exp((-rho)*(zeta[j]+1)) );
//===============================================================================


//===============================================================================
// Define seasaonality component
    season_q[j] = append_row( season_raw[j], -sum(season_raw[j]) );
    season_n[j,1:N] = season_q[j,q[1:N]];
//===============================================================================
    }
}
model {

  // priors 
  mu_u1 ~ normal(0, 1); //student_t(3, 0, 1);//
//  mu_alpha ~ normal(0, 1);//student_t(3, 0, 1);//
  mu_season ~ normal(0, 1);
  sd_season ~ normal(0, 1);//exponential(1); //
  beta ~ normal(0.5, 0.5);
  sd_y ~ normal(0, 1); //sd_y ~ exponential(1); //sd_y ~ student_t(3, 0, 1); 




  // priors for non-negative, eg, sd
//rho_raw ~ normal(0.5, 0.5);//student_t(3, 0, 1);//normal(0, 10);//1); // 
alpha_ratio_raw ~ normal(0,5);// ~ normal(0.5, 0.5); 
p ~ normal(0, 5);//5);
sd_pi ~ normal(0, 1); 

  for (j in 1:J) { 
    
//    err_pi[j] ~ normal(0,1);

    season_raw[j] ~ normal(mu_season, sd_season);
    r[j,1] ~ normal(season_n[j,1] + mu_u1, sd_y);
    r[j,2:N] ~ normal(season_n[j,2:N] + alpha_ratio*realLT[j,1:(N-1)] + beta*r[j,1:(N-1)], sd_y); 
    }

}
generated quantities {
}

