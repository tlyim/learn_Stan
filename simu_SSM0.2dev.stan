data {
  int<lower=0> J; // number of firms
  int<lower=0> N; // number of observations (for different years)
  int<lower=0> Q; // number of seasons, eg, 4 for quarters
  int<lower=0> S; // number of coefficents for the column vectors in covariate matrix T  
  int<lower=0> I; // number of coefficents for the column vectors in covariate matrix Z (which can equal X)
  int<lower=0> K; // number of coefficents for the column vectors in covariate matrix X
  int<lower=0> H; // number of coefficents for the column vectors in covariate matrix G
  int<lower=0> L; // length of simplex vector representing reversal lasting L periods
  matrix[N,S] T[J]; // covariates (capturing the dis/incentives of real earnings management)
  matrix[N,I] Z[J]; // covariates (capturing the dis/incentives of real earnings management)
  matrix[N,K] X[J]; // covariates (capturing the dis/incentives of misreporting)
  matrix[N,H] G[J]; // covariates (capturing the strength of internal governance and external monitoring)

  real<lower=0> rho; // coeff. of the LT impact of real EM
  real y_init;
  real<lower=0,upper=1>  mu_alpha; // group intercept coeff. of the AR(1) process of the underlying unbiased accounting figure y[n]
  real<lower=0,upper=1> beta; // slope coefficient of the AR(1) process of y[n]
real<lower=0> theta;
vector[S] s; // coefficients of the S covariates in matrix T
  vector[I] p; // coefficients of the K covariates in matrix Z
  vector[K] g; // coefficients of the K covariates in matrix X
  vector[H] w; // coefficients of the H covariates in matrix G
  simplex[L] d; // fractional reversal of prior-period manipluation by accruals

  vector[Q-1] mu_season;
  real<lower=0> sd_season; // sd of the hyperprior for pi
real<lower=0> sd_alpha;
  real<lower=0> sd_y; // sd of the hyperprior for alpha
  real<lower=0> sd_gamma; // sd of the hyperprior for gamma
  real<lower=0> sd_omega; // sd of the hyperprior for omega
  real mu_base; // sd of the hyperprior for base
  real<lower=0> sd_base; // sd of the hyperprior for base
}
transformed data {
}
parameters {}
model {}
generated quantities {
  vector[N] r[J]; // underlying unbiased figures (eg, Gross Profit, scaled by Rev)
  vector[N] y[J]; // underlying unbiased figures (eg, Gross Profit, scaled by Rev)
  vector[N] u[J]; // underlying unbiased figures (eg, Gross Profit, scaled by Rev)
  vector[N] Real[J]; // LT impact of real EM on alpha[j,n]
vector[N] sigma[J]; // fraction of real EM that is shifting from next period's sales
  // sigma = 0 means  the real EM is either 
  //    (i) purely RnD-based or 
  //    (ii) sales-based but not by shifting next period's sales forward

  real<lower=0> base[J];  // ,upper=1
  vector[N] alpha[J]; // 
  vector[N] b[J]; // bias effort driven by the temptation to misreport //<lower=-1,upper=1>
  vector[N] R[J];  // potential room of manipulation constrained by governance mechanisms //<lower=0>
  vector[N] m[J]; // misreporting extent in the reported figures 
  vector[N] D[J];

  vector[N] season_n[J];
  vector[Q] season_q[J];
  vector[Q-1] season_raw[J];

  // Set the quarter index
  int<lower=1,upper=Q> q[N]; // index for seasons, eg, q[n] = 1 if n belongs to Q1
  int shift = poisson_rng(1) % Q;
  for (n in 1:N) q[n] = 1 + ((n + shift) % Q); //collapse 1:N to 1:Q


// Define seasaonality component
  for (j in 1:J) {
    for (ss in 1:(Q-1)) season_raw[j,ss] = mu_season[ss] + sd_season*normal_rng(0,1);
    season_q[j] = append_row( season_raw[j], -sum(season_raw[j]) );
    }
    
  season_n[1:J,1:N] = season_q[1:J,q[1:N]];


  for (j in 1:J) {
    vector[N] tau[J]; // temptation to misreport
    vector[N] chi[J]; // 
vector[N] zeta[J]; // temptation to manage current-period real earnings upward
    vector[N] phi[J]; // 

//===============================================================================
// Define real EM component
zeta[j] = T[j]*s;// temptation to manage current-period real earnings upward
sigma[j] = ( log1p_exp(rho*zeta[j]) - log1p_exp(rho*(zeta[j]-1)) )/rho;
//===============================================================================
// Define real EM component
    phi[j] = Z[j]*p;// temptation to manage current-period real earnings upward
    Real[j] = ( log1p_exp(rho*phi[j]) - log1p_exp(rho*(phi[j]-1)) )/rho;
//===============================================================================


// fraction of Real[j] that is sales-based, rather than RnD-based,
//   where T indicates whether j is a retailer and has RnD spending or not reported in last period
//sigma[j] = inv_logit(T[j]*s);// + err_sigma[j]); 
//theta[j] = beta + sigma[j];


//===============================================================================
// Define real EM's LT impact on alpha
    alpha[j,1] = mu_alpha;  
    u[j,1] = season_n[j,1] + y_init + sd_y*normal_rng(0,1);   // y should be nonnegative for Rev and nonpositive for Costs
    for (n in 2:N) {
      alpha[j,n] = alpha[j,n-1] - theta*Real[j,n-1];  
      u[j,n] = season_n[j,n] + alpha[j,n] + beta*u[j,n-1] - sigma[j,n-1] .* Real[j,n-1] + sd_y*normal_rng(0,1);
      }

    y[j] = Real[j] + u[j]; //assumed all real EM are RnD-based; thus, no shifted sales


//Work>>> ===============================================================================
//    u[j,1] = season_n[j,1] + y_init;   // y should be nonnegative for Rev and nonpositive for Costs
//    y[j,1] = Real[j,1] + u[j,1] + sd_y*normal_rng(0,1);   // y should be nonnegative for Rev and nonpositive for Costs
//    for (n in 2:N) {
//      u[j,n] = season_n[j,n] + mu_alpha + beta*(y[j,n-1] - Real[j,n-1]);
//      y[j,n] = Real[j,n] + u[j,n] + sd_y*normal_rng(0,1); //assumed all real EM are RnD-based; thus, no shifted sales
//      //y[j,n] = Real[j,n] + u[j,n] - (theta[j,n] - beta) .* Real[j,n-1] + sd_y*normal_rng(0,1);
//      }
//Work<<< ===============================================================================


//===============================================================================
// Define accrual-based EM component
// Define the raw (m) and net (D) accrual-based EM

            //    mu_base[j] = (1)*mean( fabs(y[j]) );  // fabs() is floating-point abs()
            //    for (n in 1:N) {
            //      base[j,n] = exp( mu_base[j] + normal_rng(0, sd_base) );// sd_omega*err_omega[j]) ; // for debugging only
            //      base[j] = exp( 0.2 + 0.1*err_log_base[j] );
            //      }

    tau[j] = X[j]*g;// temptation to misreport
    chi[j] = -G[j]*w;// 

    base[j] = exp( mu_base + sd_base*normal_rng(0,1) ); 
    b[j] = 2*( log1p_exp(rho*tau[j]) - log1p_exp(rho*(tau[j]-1)) )/rho - 1;
    R[j] = ( log1p_exp(rho*chi[j]) - log1p_exp(rho*(chi[j]-1)) )/rho;

    m[j] = base[j] * b[j] .* R[j];  // extent of misreporting resulting from bias effort and 
//===============================================================================
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
//===============================================================================

    r[j] = y[j] + D[j];

    }

    
//} //end of temporary block     
}


