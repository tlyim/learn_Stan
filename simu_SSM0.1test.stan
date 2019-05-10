data {
  int<lower=0> J; // number of firms
  int<lower=0> N; // number of observations (for different years)
  int<lower=0> Q; // number of seasons, eg, 4 for quarters
  int<lower=0> I; // number of coefficents for the column vectors in covariate matrix Z (which can equal X)
  int<lower=0> K; // number of coefficents for the column vectors in covariate matrix X
  int<lower=0> H; // number of coefficents for the column vectors in covariate matrix G
  matrix[N,I] Z[J]; // covariates (capturing the dis/incentives of real earnings management)
  matrix[N,K] X[J]; // covariates (capturing the dis/incentives of misreporting)
  matrix[N,H] G[J]; // covariates (capturing the strength of internal governance and external monitoring)

real<lower=0> rho_ST; // coeff. of the ST impact of real EM
real<lower=0> rho_LT; // coeff. of the LT impact of real EM
  real y_init;
real<lower=0,upper=1>  mu_alpha; // group intercept coeff. of the AR(1) process of the underlying unbiased accounting figure y[n]
  real<lower=0,upper=1> beta; // slope coefficient of the AR(1) process of y[n]
  vector[I] p; // coefficients of the K covariates in matrix Z
  vector[K] g; // coefficients of the K covariates in matrix X
  vector[H] w; // coefficients of the H covariates in matrix G

vector[Q-1] mu_season;
real<lower=0> sd_season; // sd of the hyperprior for pi
  real<lower=0> sd_y; // sd of the hyperprior for alpha
  real<lower=0> sd_pi; // sd of the hyperprior for pi
  real<lower=0> sd_gamma; // sd of the hyperprior for gamma
  real<lower=0> sd_omega; // sd of the hyperprior for omega
  real<lower=0> sd_base; // sd of the hyperprior for base
}
transformed data {
}
parameters {}
model {}
generated quantities {
//  vector[N] alpha[J]; // LT profitability level coefficient for firm j in period n
  vector[N] y[J]; // underlying unbiased figures (eg, Gross Profit, scaled by Rev)
  vector[N] u[J]; // underlying unbiased figures (eg, Gross Profit, scaled by Rev)
  vector[N] RealLT[J]; // LT impact of real EM on alpha[j,n]


  vector[N] season_n[J];
  vector[Q] season_q[J];
  vector[Q-1] season_raw[J];

  // Set the quarter index
  int<lower=1,upper=Q> q[N]; // index for seasons, eg, q[n] = 1 if n belongs to Q1

  int shift = poisson_rng(1) % Q;
  for (n in 1:N) q[n] = 1 + ((n + shift) % Q); //collapse 1:N to 1:Q


// Define seasaonality component
  for (j in 1:J) {
    for (s in 1:(Q-1)) season_raw[j,s] = mu_season[s] + sd_season*normal_rng(0,1);
    season_q[j] = append_row( season_raw[j], -sum(season_raw[j]) );
    }
    
  season_n[1:J,1:N] = season_q[1:J,q[1:N]];



  for (j in 1:J) {
    vector[I] err_pi[J];  // error term in p
    vector[N] zeta[J]; // temptation to manage current-period real earnings upward

    for (i in 1:I) err_pi[j,i] = normal_rng(0,1);

    zeta[j] = Z[j]*(p + sd_pi*err_pi[j]);  // temptation to manage current-period real earnings upward
    RealLT[j] = ( log1p_exp((-rho_LT)*zeta[j]) - log1p_exp((-rho_LT)*(zeta[j]+1)) )/rho_LT;


    u[j,1] = season_n[j,1] + y_init + sd_y*normal_rng(0,1);   // y should be nonnegative for Rev and nonpositive for Costs
    for (n in 2:N) {
      u[j,n] = season_n[j,n] + mu_alpha*RealLT[j,n-1] + beta*u[j,n-1] + sd_y*normal_rng(0,1);
//      u[j,n] = season_n[j,n] + mu_alpha + beta*u[j,n-1] + sd_y*normal_rng(0,1);
      }


    for (n in 1:N) {
      y[j,n] = u[j,n]; 
      }

    }

    
//} //end of temporary block     
}


