// The observations are reported accounting figures, say Gross Profit, (scaled by Rev, total revenue) a vector r of length N (years).
// The covariates are a matrix capital X with K column vectors, incl a unit vector, capturing the  
//   dis/incentives of misreporting the underlying unbiased figures (y) as the observed reported firgures
data {
  int<lower=0> J; // number of firms
  int<lower=0> N; // number of observations (for different years)
  int<lower=0> Q; // number of seasons, eg, 4 for quarters
  int<lower=1,upper=Q> q[N]; // index for seasons, eg, q[n] = 1 if n belongs to Q1
  int<lower=0> I; // number of coefficents for the column vectors in covariate matrix Z (which can equal X)
  int<lower=0> K; // number of coefficents for the column vectors in covariate matrix X
  int<lower=0> H; // number of coefficents for the column vectors in covariate matrix G
  int<lower=0> L; // length of simplex vector representing reversal lasting L periods

real mu_alpha01_init;
real beta_init;
real delta_init;
real rho_ST_init;
real rho_LT_init;
  vector[N] r[J]; // reported figures (scaled by Rev, total revenue)
  matrix[N,I] Z[J]; // covariates (capturing the dis/incentives of real earnings management)
  matrix[N,K] X[J]; // covariates (capturing the dis/incentives of misreporting)
  matrix[N,H] G[J]; // covariates (capturing the strength of internal governance and external monitoring mechansims)
//  real<lower=0> base;
// forecasts  
//  int<lower=0> N_new; // number of predictions
//  matrix[N_new,K] X_new; // 
}
transformed data {
}
parameters {
  real<lower=0,upper=1> mu_alpha01; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
  real<lower=0,upper=1> beta; // slope coefficient of the AR(1) process of y[n]
//  real<lower=0,upper=1> delta; // slope coefficient of the AR(1) process of y[n]
  real<lower=0> delta; // slope coefficient of the AR(1) process of y[n]
real<lower=0> rho_ST; // coeff. of the ST impact of real EM
real<lower=0> rho_LT; // coeff. of the LT impact of real EM
  vector[I] p; // coefficients of the K covariates in matrix Z
  vector[K] g; // coefficients of the K covariates in matrix X
  vector[H] w; // coefficients of the H covariates in matrix G
  simplex[L] d; // fractional reversal of prior-period manipluation by accruals

  real<lower=0> sd_y; // sd of the underlying unbiased figures (vector y)
//  real<lower=0> sd_alpha; // sd of the underlying unbiased figures (vector y)
  real<lower=0> sd_pi; // sd of the hyperprior for pi
  real<lower=0> sd_gamma; // sd of the hyperprior for gamma
  real<lower=0> sd_omega; // sd of the hyperprior for omega
  real<lower=0> sd_base; // sd of the hyperprior for gamma
  real<lower=0,upper=1> mu_base[J]; // sd of the hyperprior for gamma
  vector[N] err_log_base[J];

  vector[I] err_pi[J];  // error term in p
//vector[N] err_alpha[J];  // error term in p
  vector[H] err_omega[J];
  vector[K] err_gamma[J];
//vector[N] err_season[J];  
//real<lower=0> sd_season; // sd of the underlying unbiased figures (vector y)
vector[Q-1] season_raw[J];
vector[Q-1] mu_season;
real<lower=0> sd_season; // sd of the hyperprior for pi

// real<lower=0,upper=1> integrity;     // integrity level of the society affecting the decision to misreport or not
//  vector[N] integrity;     // integrity level of each CEO affecting the decision to misreport or not
}
transformed parameters {
  vector[N] y[J];
  vector[N] m[J]; // misreporting extent in the reported figures 
vector[N] RealST[J]; // ST impact of real EM on y[j,n]
vector[N] alpha[J]; // underlying unbiased figure (scaled by Rev, total revenue)
  real<lower=-1,upper=1> mu_alpha; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]

// Define seasaonality component
vector[Q] season_q[J];
vector[N] season_n[J];

  for (j in 1:J) {
    season_q[j] = append_row( season_raw[j], -sum(season_raw[j]) );
    for (n in 1:N) {   
      season_n[j,n] = season_q[j,q[n]];
      }
    }

// Define (-1,1)-bounded mean alpha as a transformation of the (0,1)-bounded mean alpha01
mu_alpha = 2*mu_alpha01 - 1;


// Define the components of real and accrual-based EM 
  for (j in 1:J) {  
  vector[N] D[J];
    vector[N] tau[J]; // temptation to misreport
vector[N] chi[J]; // 
vector[N] base[J];  // ,upper=1
vector[N] zeta[J]; // temptation to manage current-period real earnings upward
vector[N] RealLT[J]; // LT impact of real EM on alpha[j,n]


// Define the sT and LT impacts of real EM 
   zeta[j] = Z[j]*(p + sd_pi*err_pi[j]);  // temptation to manage current-period real earnings upward
RealST[j] = rep_vector(rho_ST, N) ./ ( 1 + exp((-1)*zeta[j]) );
//RealST[j] = rep_vector(0.1, N) ./ ( 1 + exp((-1)*zeta[j]) );
//RealLT[j] = rep_vector(0.5*0.04, N) ./ ( 1 + exp((-1)*zeta[j]) );
RealLT[j] = rep_vector(rho_LT, N) ./ ( 1 + exp((-1)*zeta[j]) );


// Define the raw (m) and net (D) accrual-based EM
    tau[j] = X[j]*(g + sd_gamma*err_gamma[j]);  // temptation to misreport
    chi[j] = G[j]*(w + sd_omega*err_omega[j]);
base[j] = exp( mu_base[j] + sd_base*err_log_base[j] );

//    vector[N] b[J]; // bias effort driven by the temptation to misreport //<lower=-1,upper=1>
//    vector[N] P[J];  // potential room of manipulation constrained by governance mechanisms //<lower=0>
 //    b[j] = (exp(tau[j]) - 1) ./ (exp(tau[j]) + 1);
//    P[j] = inv_logit( (-1)*G[j]*(w + sd_omega*err_omega[j]) );

//    m[j] = b[j] .* P[j] .* rep_vector(base[j], N);  // extent of misreporting resulting from bias effort 
    m[j] = ( base[j] .* (exp(tau[j]) - 1) ) ./ ( 1 + exp(tau[j]) + exp(chi[j]) + exp(tau[j]+chi[j]) );

    
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
    for (l in 1:L) {
//      D[j,(L+1):N] += (-1)* d[l] * m[j,(L+1-l):(N-l)];
      D[j,(L+1):N] = D[j,(L+1):N] + (-1)* d[l] * m[j,(L+1-l):(N-l)];
      }


    y[j] = r[j] - D[j];
// test this again after RealST is included:    
    alpha[j,1] = mu_alpha - RealLT[j,1];// + sd_alpha * err_alpha[j,1];
//    alpha[j,2:N] = alpha[j,1:(N-1)];// - RealLT[j,2:N];// + sd_alpha * err_alpha[j,2:N];
    for (n in 2:N) { 
      alpha[j,n] = alpha[j,n-1] - RealLT[j,n];// + sd_alpha * err_alpha[j,n];
      }

    
    }

//  real ilogodds;     // log odds of the society integrity level
//  ilogodds = logit(integrity);   
}
model {

// priors 
//	sd_y ~ cauchy(0, 5);//normal(0, 1);//
//sd_pi ~ cauchy(0, 5);//normal(0, 1);//
//	sd_gamma ~ cauchy(0, 5);//normal(0, 1);//
//	sd_omega ~ cauchy(0, 5);//normal(0, 1);//
//  sd_base ~ cauchy(0, 5);//normal(0, 1);//
//delta ~ normal(0, 10);//cauchy(0, 5);//
  	sd_y ~ normal(0, 10);//cauchy(0, 5);//
  sd_base ~ normal(0, 10);//cauchy(0, 5);//
// 	sd_gamma ~ normal(0, 10);//cauchy(0, 5);//
// 	sd_omega ~ normal(0, 10);//cauchy(0, 5);//
// sd_pi ~ normal(0, 10);//cauchy(0, 5);//
  

delta ~ normal(0, 5);//2.5);//cauchy(0, 5);//
rho_ST ~ normal(0, 0.5);//0.25);// 0.125);// uniform(0,0.2);
//rho_ST ~ exponential(0.1);// uniform(0,0.2);
//rho_LT ~ exponential(20);// uniform(0,0.2);
rho_LT ~ normal(0, 0.5);//0.25);//0.125);// uniform(0,0.2);


sd_pi ~ normal(0, 0.125);//cauchy(0, 5);//
	sd_gamma ~ normal(0, 0.125);//cauchy(0, 5);//
	sd_omega ~ normal(0, 0.125);//cauchy(0, 5);//
sd_season ~ normal(0, 0.125);//cauchy(0, 5);//
//  mu_alpha01 ~ beta(3*mu_alpha01_init, 3*(1-mu_alpha01_init));
//  beta ~ beta(5.5*beta_init, 5.5*(1-beta_init));
  mu_season ~ normal(0, 10); //  10); //0.2); //  
  p ~ normal(0, 10); //  10); //0.2); //  
  g ~ normal(0, 10);//0); //0.2); // g_init
  w ~ normal(0, 10);//0); //0.2); //1);  // w_init
  
//  d ~ dirichlet(rep_vector(0.1, L));   // =1.0 means the dist is uniform over the possible simplices;<1.0 toward corners 



  for (j in 1:J) { 

//    err_alpha[j] ~ normal(0, 1);
//    err_season[j] ~ normal(0, 1);
    err_pi[j] ~ normal(0, 1);
    err_gamma[j] ~ normal(0, 1);
    err_omega[j] ~ normal(0, 1);
    err_log_base[j] ~ normal(0, 1);

//=======================================================================================
// Two options: 
//(i) simply let y[j,1] be an unknown parameter to be best estimated to fit the data; a prior can be specified
//y[j,1] ~ normal(0, 1);
//
//(ii) model y[j,1] explicitly, perhaps with a measure error term to capture the unknown - delta*RealST[j,1:(N-1)] + RealST[j,2:N]
//y[j,1] ~ normal( season_n[j,1] + sd_season * err_season[j,1] + mu_alpha/(1-beta), sd_y );
//=======================================================================================

//y[j,2:N] = alpha[j,1:(N - 1)] + beta*(y[j,1:(N - 1)] - RealST[j,1:(N - 1)]) + RealST[j,2:N]; 
//      y[j,2:N] ~ normal(mu_alpha + beta *( y[j,1:(N-1)] - RealST[j,1:(N - 1)] ) + RealST[j,2:N], sd_y);
season_raw[j] ~ normal(mu_season, sd_season);
      y[j,2:N] ~ normal( season_n[j,2:N] + //sd_season * err_season[j,2:N]
                          alpha[j,1:(N-1)] + beta*y[j,1:(N-1)] +
                          delta*RealST[j,1:(N-1)] + RealST[j,2:N], sd_y );
//      y[j,2:N] ~ normal(alpha[j,1:N-1] + beta * y[j,1:N-1], sd_alpha);

  }
//  M ~ bernoulli_logit(ilogodds); 
}
generated quantities {

//  vector[N_new] r_new;
//  for (n in 1:N_new)
//    r_new[n] = normal_rng(Z_new[n] * b, sd_r);
}

