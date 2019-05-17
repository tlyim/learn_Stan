data {
  int<lower=0> J; // number of firms
  int<lower=0> N; // number of observations (for different years)
  int<lower=0> Q; // number of seasons, eg, 4 for quarters
  int<lower=1,upper=Q> q[N]; // index for seasons, eg, q[n] = 1 if n belongs to Q1
  int<lower=0> S; // number of coefficents for the column vectors in covariate matrix T
  int<lower=0> I; // number of coefficents for the column vectors in covariate matrix Z (which can equal X)
  int<lower=0> K; // number of coefficents for the column vectors in covariate matrix X
  int<lower=0> H; // number of coefficents for the column vectors in covariate matrix G
  int<lower=0> L; // length of simplex vector representing reversal lasting L periods

real rho; // coeff. of the LT impact of real EM
  vector[N] r[J]; // reported figures (scaled by Rev, total revenue)
  matrix[N,S] T[J]; // covariates (capturing the dis/incentives of real earnings management)
  matrix[N,I] Z[J]; // covariates (capturing the dis/incentives of real earnings management)
  matrix[N,K] X[J]; // covariates (capturing the dis/incentives of misreporting)
  matrix[N,H] G[J]; // covariates (capturing the strength of internal governance and external monitoring mechansims)
}
transformed data {
//  real rho; // coeff. of the LT impact of real EM
//rho = 1;

}
parameters {
  real mu_u1;
  real mu_alpha; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
//real alpha_ratio_raw; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
  real<lower=0,upper=1> beta; // slope coefficient of the AR(1) process of y[n]

//  real<lower=0,upper=1> eta; // slope coefficient of the AR(1) process of y[n]
  real<lower=0> sd_y; // sd of the underlying unbiased figures (vector y)
  vector[Q-1] season_raw[J];
  vector[Q-1] mu_season;   // mu of the hyperprior for season_raw
  real<lower=0> sd_season; // sd of the hyperprior for season_raw

//real<lower=0,upper=1> rho_raw; // coeff. of the LT impact of real EM
//real rho_raw; // coeff. of the LT impact of real EM
//vector[I-1] p_raw; // coefficients of the K covariates in matrix Z
  vector[I] p; // coefficients of the K covariates in matrix Z
//real<lower=0> sd_pi; // sd of the hyperprior for pi
//vector[N] err_pi[J];  // error term in p
  vector[S] s; // coefficients of the S covariates in matrix T
}
transformed parameters {
//vector[N] u[J];
//vector[N] err[J];  // error term in p
  vector[Q] season_q[J];
  vector[N] season_n[J];
vector[N] Real[J]; // LT impact of real EM on alpha[j,n]

vector<lower=0,upper=1>[N] sigma[J]; // slope coefficient of the AR(1) process of y[n]
vector[N] theta[J]; // slope coefficient of the AR(1) process of y[n]
//real mu_alpha; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]

  for (j in 1:J) {
    vector[N] zeta[J]; // temptation to manage current-period real earnings upward
//===============================================================================
// Define real EM component
    zeta[j] = Z[j]*p;// + sd_pi*err_pi[j];  // temptation to manage current-period real earnings upward
//    RealLT[j] = log( log1p_exp((-rho)*zeta[j]) - log1p_exp((-rho)*(zeta[j]+1)) ) - log(rho);
    Real[j] = ( log1p_exp(rho*zeta[j]) - log1p_exp(rho*(zeta[j]-1)) )/rho;
//    Real[j] = (0.1)*inv_logit(zeta[j]);

//    realLT[j] = ( log1p_exp((-1)*zeta_rho[j]) - log1p_exp((-1)*(zeta_rho[j]+rho)) ) ;
//    realLT[j] = ( log1p_exp((-rho)*zeta[j]) - log1p_exp((-rho)*(zeta[j]+1)) ) ;
//===============================================================================


//===============================================================================
// Define seasaonality component
    season_q[j] = append_row( season_raw[j], -sum(season_raw[j]) );
    season_n[j,1:N] = season_q[j,q[1:N]];
//===============================================================================


// fraction of Real[j] that is sales-based, rather than RnD-based,
//   where T indicates whether j is a retailer and has RnD spending or not reported in last period
sigma[j] = inv_logit(T[j]*s);// + err_sigma[j]); 

theta[j] = beta + sigma[j];
//theta[j] = beta + sigma[j] + eta*(1-sigma[j]);

    }
}
model {

  // priors 
  mu_u1 ~ normal(0, 1); //student_t(3, 0, 1);//
  mu_alpha ~ normal(0, 1);//student_t(3, 0, 1);//
  mu_season ~ normal(0, 1);
  sd_season ~ normal(0, 1);//exponential(1); //
  beta ~ normal(0.5, 0.5);
  sd_y ~ normal(0, 1); //sd_y ~ exponential(1); //sd_y ~ student_t(3, 0, 1); 

  // priors for non-negative, eg, sd

p ~ student_t(4, 0, 1);//normal(0, 1);//5);
s ~ normal(0, 1);//5);//student_t(4, 0, 1);//
//sd_pi ~ normal(0, 1); 
//sigma ~ beta(0.5,0.5);
//shape ~ normal(0, 5);//student_t(4, 0, 1);// normal(0, 1); //cauchy(0,5); //half-normal
//eta ~ normal(0.5, 0.5);


  for (j in 1:J) { 
    
  for (n in 1:N)   

    season_raw[j] ~ normal(mu_season, sd_season);
    r[j,1] ~ normal(Real[j,1] + season_n[j,1] + mu_u1, sd_y);
//    r[j,2:N] ~ normal(Real[j,2:N] + season_n[j,2:N] + mu_alpha + beta*(r[j,1:(N-1)] - Real[j,1:(N-1)]), sd_y); 
    r[j,2:N] ~ normal(Real[j,2:N] + season_n[j,2:N] - theta[j,2:N] .* Real[j,1:(N-1)] +
                      mu_alpha + beta*r[j,1:(N-1)], 
                      sd_y); 
//    r[j,2:N] ~ normal(Real[j,2:N] + season_n[j,2:N] - sigma*Real[j,1:(N-1)] +
//                      mu_alpha - eta*(1-sigma)*Real[j,1:(N-1)] +
//                      beta*(r[j,1:(N-1)] - Real[j,1:(N-1)]), 
//                      sd_y); 

//Given Real[j,n-1] of last period 
//  and that sigma in [0,1] is the fraction driven by speeding up sales to last period from this period's, 
//  with (1 - sigma) being the fraction of real EM driven by cutting RnD spending,
// the impacts of last period's real EM in this period are:
//   (i) the kernel earnings of last period forming the basis for determining the persistent earnings
//       should be (r[j,n-1] - Real[j,n-1]) only;
//   (ii) moreover, the kernel earnings of this year would be sigma*Real[j,n-1]) less due to 
//       the loss of sales that had been shifted forward to last period (would be zero if signma = 0);
//   (iii) furthermore, the RnD-driven fraction of real EM would reduce the LT profitability 
//       through reducing mu_alpha by an amount eta*(1-sigma)*Real[j,n-1], where eta in (0,1).
// sigma is assumed to be a r.v. following beta(0.5,0.5), which has a U-shaped density, which
//   represents the view that firms tend to use either channel instead of a mix of both 
// In general, a hyperprior can be put on sigma to let firm characteristics drive the choice, eg,
//   retailers are unlikely to have a lot of RnD spending to cut. 
    }

}
generated quantities {
}

