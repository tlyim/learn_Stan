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
//vector<lower=-1,upper=1>[N] b[J]; // bias effort driven by the temptation to misreport //<lower=-1,upper=1>
//for (j in 1:J) {
//  for (n in 1:N) {
//b[j,n] = 2*inv_logit(normal_rng(0,1)) - 1;
//    }
//  }
//
//vector[N] tau[J]; // temptation to misreport
//vector[N] chi[J]; // 
//  for (j in 1:J) {
//    tau[j] = X[j]*(g);// + sd_gamma*err_gamma[j]);  // temptation to misreport
//    chi[j] = G[j]*(w);// + sd_omega*err_omega[j]);
//    }
}
parameters {
  real mu_u1;
  real mu_alpha; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
  real<lower=0,upper=1> beta; // slope coefficient of the AR(1) process of y[n]

  real<lower=0> sd_y; // sd of the underlying unbiased figures (vector y)
  vector[Q-1] season_raw[J];
  vector[Q-1] mu_season;   // mu of the hyperprior for season_raw
  real<lower=0> sd_season; // sd of the hyperprior for season_raw


//  vector[S] s; // coefficients of the S covariates in matrix T
  vector[I] p; // coefficients of the K covariates in matrix Z
vector[K] g; // coefficients of the K covariates in matrix X
vector[H] w; // coefficients of the H covariates in matrix G
simplex[L] d; // fractional reversal of prior-period manipluation by accruals
//
//
real<lower=0> sd_base; // sd of the hyperprior for gamma
real mu_base;  
real err_log_base[J];

}
transformed parameters {
  vector[Q] season_q[J];
  vector[N] season_n[J];
  vector[N] Real[J]; // LT impact of real EM on alpha[j,n]

vector<lower=-1,upper=1>[N] b[J]; // bias effort driven by the temptation to misreport //<lower=-1,upper=1>
vector<lower=0,upper=1>[N] P[J];  // potential room of manipulation constrained by governance mechanisms //<lower=0>
vector[N] m[J]; // misreporting extent in the reported figures 
vector[N] D[J];
vector[N] u[J]; // unmanaged earnings (if shock removed, the remaining is the kernel earnings)


//vector<lower=0,upper=1>[N] sigma[J]; // slope coefficient of the AR(1) process of y[n]
//vector[N] theta[J]; // slope coefficient of the AR(1) process of y[n]


  for (j in 1:J) {

vector[N] tau[J]; // temptation to misreport
vector[N] chi[J]; // 
real base[J];  // ,upper=1

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
    //sigma[j] = inv_logit(T[j]*s);// + err_sigma[j]); 
    //theta[j] = beta + sigma[j];
        //theta[j] = beta + sigma[j] + eta*(1-sigma[j]);


//===============================================================================
// Define accrual-based EM component
// Define the raw (m) and net (D) accrual-based EM

    tau[j] = X[j]*(g);// + sd_gamma*err_gamma[j]);  // temptation to misreport
    chi[j] = -G[j]*(w);// + sd_omega*err_omega[j]);

//base[j] = exp( mu_base[j] + sd_base*err_log_base[j] );
base[j] = exp( mu_base + sd_base*err_log_base[j] );//rep_vector(0.2, N);//

b[j] = 2*( log1p_exp(rho*tau[j]) - log1p_exp(rho*(tau[j]-1)) )/rho - 1;
//    b[j] = (exp(tau[j]) - 1) ./ (exp(tau[j]) + 1);
//    P[j] = inv_logit( (-1)*G[j]*(w) );  // + sd_omega*err_omega[j]
//    P[j] = inv_logit( (-1)*G[j]*(w) );  // + sd_omega*err_omega[j]
P[j] = ( log1p_exp(rho*chi[j]) - log1p_exp(rho*(chi[j]-1)) )/rho;

    m[j] = base[j] * b[j] .* P[j];  // extent of misreporting resulting from bias effort 
//    m[j] = ( base[j] .* (exp(tau[j]) - 1) ) ./ ( 1 + exp(tau[j]) + exp(chi[j]) + exp(tau[j]+chi[j]) );
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

    u[j] = r[j] - D[j] - Real[j];

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

sd_base ~ normal(0, 1);
mu_base ~ normal(0, 1);
      //s ~ normal(0, 1);//5);//student_t(4, 0, 1);//
      //sigma ~ beta(0.5,0.5);
p ~ normal(0, 1);//student_t(4, 0, 1);//5);
g ~ normal(0, 1);//student_t(4, 0, 1);//student_t(3, 0, 5);//normal(0, 1);//10);//0); //0.2); // g_init
w ~ normal(0, 1);//student_t(4, 0, 1);//student_t(3, 0, 5);//normal(0, 5);//10);//0); //0.2); //1);  // w_init
//  d ~ dirichlet(rep_vector(0.1, L));   // =1.0 means the dist is uniform over the possible simplices;<1.0 toward corners 


  for (j in 1:J) { 

// ??? identifiable ??    
//mu_base[j] ~ normal(0, 1);

err_log_base[j] ~ normal(0, 1);
    season_raw[j] ~ normal(mu_season, sd_season);

    u[j,1] ~ normal(season_n[j,1] + mu_u1, sd_y);
    u[j,2:N] ~ normal(season_n[j,2:N] + mu_alpha + beta*u[j,1:(N-1)], sd_y); 


    
//Also Work>>> ===============================================================================
//    r[j,1] ~ normal(0*D[j,1] + Real[j,1] + season_n[j,1] + mu_u1, sd_y);
//    r[j,2:N] ~ normal(0*D[j,2:N] + Real[j,2:N] + season_n[j,2:N] + //- theta[j,2:N] .* Real[j,1:(N-1)] +
//                      mu_alpha + beta*( r[j,1:(N-1)] - 0*D[j,1:(N-1)] - Real[j,1:(N-1)] ), 
//                      sd_y); 
//Also Work>>> ===============================================================================

    
//Work>>> ===============================================================================
//    r[j,1] ~ normal(0*D[j,1] + Real[j,1] + season_n[j,1] + mu_u1, sd_y);
          //    r[j,2:N] ~ normal(Real[j,2:N] + season_n[j,2:N] + mu_alpha + beta*(r[j,1:(N-1)] - Real[j,1:(N-1)]), sd_y); 
//    r[j,2:N] ~ normal(0*D[j,2:N] + Real[j,2:N] + season_n[j,2:N] + //- theta[j,2:N] .* Real[j,1:(N-1)] +
//                                    mu_alpha + beta*(r[j,1:(N-1)] - 0*D[j,1:(N-1)] - Real[j,1:(N-1)]), 
//                                    sd_y); 
                //    r[j,2:N] ~ normal(Real[j,2:N] + season_n[j,2:N] - theta[j,2:N] .* Real[j,1:(N-1)] +
                //                      mu_alpha + beta*r[j,1:(N-1)], 
                //                      sd_y); 
                //    r[j,2:N] ~ normal(Real[j,2:N] + season_n[j,2:N] - sigma*Real[j,1:(N-1)] +
                //                      mu_alpha - eta*(1-sigma)*Real[j,1:(N-1)] +
                //                      beta*(r[j,1:(N-1)] - Real[j,1:(N-1)]), 
                //                      sd_y); 
//Work<<< ===============================================================================


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

