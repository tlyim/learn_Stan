data {
  int<lower=0> J; // number of firms
  int<lower=0> N; // number of observations (for different years)
  int<lower=0> I; // number of coefficents for the column vectors in covariate matrix Z (which can equal X)
  int<lower=0> K; // number of coefficents for the column vectors in covariate matrix X
  int<lower=0> H; // number of coefficents for the column vectors in covariate matrix G
  matrix[N,I] Z[J]; // covariates (capturing the dis/incentives of real earnings management)
  matrix[N,K] X[J]; // covariates (capturing the dis/incentives of misreporting)
  matrix[N,H] G[J]; // covariates (capturing the strength of internal governance and external monitoring)

real<lower=0> rho_ST; // coeff. of the ST impact of real EM
real<lower=0> rho_LT; // coeff. of the LT impact of real EM
  real y_init;
real mu_alpha; // group intercept coeff. of the AR(1) process of the underlying unbiased accounting figure y[n]
  real<lower=0,upper=1> beta; // slope coefficient of the AR(1) process of y[n]
real<lower=0,upper=1> delta; // coeff. of the ST impact of real EM
  vector[I] p; // coefficients of the K covariates in matrix Z
  vector[K] g; // coefficients of the K covariates in matrix X
  vector[H] w; // coefficients of the H covariates in matrix G

  real<lower=0> sd_y; // sd of the hyperprior for alpha
  real<lower=0> sd_pi; // sd of the hyperprior for pi
  real<lower=0> sd_gamma; // sd of the hyperprior for gamma
  real<lower=0> sd_omega; // sd of the hyperprior for omega
  real<lower=0> sd_base; // sd of the hyperprior for base
//  real<lower=0> sd_y; // sd of the underlying unbiased figure (vector y)
}
transformed data {

}
parameters {}
model {}
generated quantities {
vector[N] RealST[J]; // ST impact of real EM on y[j,n]
vector[N] RealLT[J]; // LT impact of real EM on alpha[j,n]
  vector[N] alpha[J]; // LT profitability level coefficient for firm j in period n
  vector[N] y[J]; // underlying unbiased figures (eg, Gross Profit, scaled by Rev)
  vector[N] m[J]; // misreporting extent in the reported figures 
//  vector[N] b[J]; //<lower=-1,upper=1> bias effort driven by the temptation to misreport
//  vector[N] bR[J]; //<lower=-1,upper=1> bias effort driven by the temptation to misreport
//  vector[N] R[J];  //<lower=0> potential room of manipulation constrained by governance mechanisms
  vector<lower=0>[N] base[J];  // ,upper=1
  real mu_base[J];  // ,upper=1

  for (j in 1:J) {
    vector[K] err_gamma[J];
    vector[H] err_omega[J];
    vector[I] err_pi[J];  // error term in p
    vector[N] err_y [J];  // error term in p
//    vector[N] RealST[J]; // ST impact of real EM on y[j,n]
//    vector[N] RealLT[J]; // LT impact of real EM on alpha[j,n]
    vector[N] zeta[J]; // temptation to manage current-period real earnings upward
    vector[N] tau[J]; // temptation to misreport
    vector[N] chi[J]; // 

    for (i in 1:I) {
      err_pi[j,i] = normal_rng(0, 1);
      }
    zeta[j] = Z[j]*(p + sd_pi*err_pi[j]);  // temptation to manage current-period real earnings upward

RealST[j] = rep_vector(rho_ST, N) ./ ( 1 + exp((-1)*zeta[j]) );
RealLT[j] = rep_vector(rho_LT, N) ./ ( 1 + exp((-1)*zeta[j]) );

// Perhaps best to model y as always positive, with IsCost = 1, 0 to indicate Cost or Rev item 

    err_y[j,1] = normal_rng(0, 1);
//    alpha[j,1] = mu_alpha + sd_y * err_alpha[j,1] - RealLT[j,1];   // y should be nonnegative for Rev and nonpositive for Costs
    alpha[j,1] = mu_alpha - RealLT[j,1];   // y should be nonnegative for Rev and nonpositive for Costs
    y[j,1] = y_init + RealST[j,1] + sd_y * err_y[j,1];   // y should be nonnegative for Rev and nonpositive for Costs
    for (n in 2:N) {
      err_y[j,n] = normal_rng(0, 1);
 //     alpha[j,n] = alpha[j,n-1] + sd_y * err_alpha[j,n] - RealLT[j,n];
      alpha[j,n] = alpha[j,n-1] - RealLT[j,n];
      y[j,n] = alpha[j,n-1] + beta*y[j,n-1] - delta*RealST[j,n-1] + RealST[j,n] + sd_y * err_y[j,n]; //- delta*RealST[j,n-1];
      }
// Initialize alpha[j,1] with group alpha (ie, mu_alpha)


// Initialize at a level three times of the LT stationary level
    for (n in 2:N) {
      }
// Perhaps best to model y as always positive, with IsCost = 1, 0 to indicate Cost or Rev item 


    mu_base[j] = (1)*mean( fabs(y[j]) ); // abs() has a kink at zero
    for (n in 1:N) {
      base[j,n] = exp( mu_base[j] + normal_rng(0, sd_base) );// + sd_omega*err_omega[j]) ; // for debugging only
      }


    for (k in 1:K) {
      err_gamma[j,k] = normal_rng(0, 1);
      }
    for (h in 1:H) {
      err_omega[j,h] = normal_rng(0, 1);
      }

    tau[j] = X[j]*(g + sd_gamma*err_gamma[j]);  // temptation to misreport
    chi[j] = G[j]*(w + sd_omega*err_omega[j]);

//    b[j] = (exp(tau[j]) - 1) ./ (exp(tau[j]) + 1);
//    R[j] = rep_vector(1, N) ./ ( 1 + exp(chi[j]) );
//    bR[j] = (exp(tau[j]) - 1) ./ ( 1 + exp(tau[j]) + exp(chi[j]) + exp(tau[j]+chi[j]) );

//    m[j] = b[j] .* R[j] .* rep_vector(base[j], N);  // extent of misreporting resulting from bias effort and 
//    m[j] = rep_vector(base[j], N) .* bR[j];  // extent of misreporting resulting from bias effort and 
    m[j] = ( base[j] .* (exp(tau[j]) - 1) ) ./ ( 1 + exp(tau[j]) + exp(chi[j]) + exp(tau[j]+chi[j]) );
    }
}


