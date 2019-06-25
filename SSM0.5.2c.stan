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
//real mu_base;
  real rho; // coeff. of the LT impact of real EM
  vector[N] r[J]; // reported figures (scaled by Rev, total revenue)
  matrix[N,S] T[J]; // covariates (capturing the dis/incentives of real earnings management)
  matrix[N,I] Z[J]; // covariates (capturing the dis/incentives of real earnings management)
  matrix[N,K] X[J]; // covariates (capturing the dis/incentives of misreporting)
  matrix[N,H] G[J]; // covariates (capturing the strength of internal governance and external monitoring mechansims)
}
transformed data {
  real p3 = (-10)*0.7; //<--------------------------------
//!!!
  int<lower=0> I_cor = I;//2; // number of correlated coefficents in p for covariate matrix Z
  int<lower=0> K_cor = K; // number of correlated coefficents in g for covariate matrix X
  int<lower=0> H_cor = H; // number of correlated coefficents in w for covariate matrix G
}
parameters {
  real<lower=0> sd_y; // sd of the underlying unbiased figures (vector y)
vector[Q-1] season_raw[J];
vector[Q-1] mu_season;   // mu of the hyperprior for season_raw
real<lower=0> sd_season; // sd of the hyperprior for season_raw

  real mu_u1;
real<lower=0> theta; //reduction in alpha[j,n] from the initial level (mu_alpha) due to real EM's LT impact 
//===============================================================================
//!!!  real<lower=0,upper=1> mu_alpha; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
//!!!  real<lower=0,upper=1> beta; // slope coefficient of the AR(1) process of y[n]
    vector[2] ab_mu; //<lower=0,upper=1> vector of means for correlated coeffs in w
    vector<lower=0>[2] ab_sd; // vector of sd for correlated coeffs in w
    cholesky_factor_corr[2] ab_L; // Cholesky factor of correlation matrix for correlated coeffs in w
    vector[2] ab_err; // primitive vector of correlated coeffs in w
//===============================================================================

vector[I-1] p_raw;
vector[K-1] g_raw;
vector[H-1] w_raw;

//===============================================================================
//  vector[K] g; // coefficients of the K covariates in matrix X
            //vector[K - K_cor] g0; // vector of means for correlated coeffs in w
//!!! Trying to require >= 0 to get ride of an unreasonable second mode             
vector<lower=0>[3] pgw_mu; //<lower=0> vector of means for correlated coeffs in w
    vector<lower=0>[3] pgw_sd; // vector of sd for correlated coeffs in w
    cholesky_factor_corr[3] pgw_L; // Cholesky factor of correlation matrix for correlated coeffs in w
    vector[3] pgw_err; // primitive vector of correlated coeffs in w
/*
real g1;
real g2;
vector[2] gw_mu;
vector<lower=0>[2] gw_sd; // vector of sd for correlated coeffs in w
cholesky_factor_corr[2] gw_L; // Cholesky factor of correlation matrix for correlated coeffs in w
vector[2] gw_err; // primitive vector of correlated coeffs in w
real w2;
real w3;
//===============================================================================
//  vector[K] g; // coefficients of the K covariates in matrix X
            //vector[K - K_cor] g0; // vector of means for correlated coeffs in w
//            vector[K_cor] g_mu; // vector of means for correlated coeffs in w
            vector<lower=0>[K_cor] g_sd; // vector of sd for correlated coeffs in w
            cholesky_factor_corr[K_cor] g_L; // Cholesky factor of correlation matrix for correlated coeffs in w
            vector[K_cor] g_err; // primitive vector of correlated coeffs in w
//===============================================================================
//  vector[H] w; // coefficients of the H covariates in matrix G
            //vector[H - H_cor] w0; // vector of means for correlated coeffs in w
//            vector[H_cor] w_mu; // vector of means for correlated coeffs in w
            vector<lower=0>[H_cor] w_sd; // vector of sd for correlated coeffs in w
            cholesky_factor_corr[H_cor] w_L; // Cholesky factor of correlation matrix for correlated coeffs in w
            vector[H_cor] w_err; // primitive vector of correlated coeffs in w
*/
//===============================================================================
  simplex[L] d; // fractional reversal of prior-period manipluation by accruals
//
//!!!  real mu_base;  
//  real<lower=0> sd_base; // sd of the hyperprior for gamma
//  real err_log_base[J];
}
transformed parameters {
//vector[Q] season_q[J];
//vector[N] season_n[J];
//!!!    vector[N] u[J]; // unmanaged earnings (if shock removed, the remaining is the kernel earnings)
//  vector[N] Real[J]; // LT impact of real EM on alpha[j,n]
//  vector[N] m[J]; // misreporting extent in the reported figures 


//vector[N] sigma[J]; // fraction of real EM that is shifting from next period's sales
  // sigma = 0 means  the real EM is either 
  //    (i) purely RnD-based or 
  //    (ii) sales-based but not by shifting next period's sales forward


//===============================================================================
  real mu_alpha; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
  real beta; // slope coefficient of the AR(1) process of y[n]
//===============================================================================
  vector[I] p; // coefficients of the H covariates in matrix G
  vector[K] g; // coefficients of the K covariates in matrix X
  vector[H] w; // coefficients of the H covariates in matrix G


//===============================================================================
{
  vector[2] ab; // primitive vector of correlated coeffs in w
  ab = ab_mu + ab_sd .* (ab_L * ab_err);
    mu_alpha = ab[1];//inv_logit(ab[1]);
    beta = ab[2];
}    
//===============================================================================
{
  vector[3] pgw; // modeling g and w together to capture their correlations
  pgw = pgw_mu + pgw_sd .* (pgw_L * pgw_err);
    p[1] = pgw[1];
    g[1] = pgw[2];
    w[1] = pgw[3];
    p[2:I] = p_raw[1:I-1];
    g[2:K] = g_raw[1:K-1];
    w[2:H] = w_raw[1:H-1];
}    
//===============================================================================


}
model {

vector[N] err_y[J];

  for (j in 1:J) {
vector[Q] season_q[J];
vector[N] season_n[J];
vector[N] Real[J]; // LT impact of real EM on alpha[j,n]
vector[N] m[J]; // misreporting extent in the reported figures 

    vector[N] alpha[J];
//real base[J];  // ,upper=1
vector[N] zeta[J]; // temptation to manage current-period real earnings upward
    vector[N] tau[J]; // temptation to misreport
    vector[N] chi[J]; // 
    vector[N] phi[J]; // 

    vector[N] b[J]; //<lower=-1,upper=1> bias effort driven by the temptation to misreport //<lower=-1,upper=1>
    vector[N] R[J];  //<lower=0,upper=1> potential room of manipulation constrained by governance mechanisms //<lower=0>
    vector[N] D[J];


//===============================================================================
// Define accrual-based EM component
// Define the raw (m) and net (D) accrual-based EM

    tau[j] = X[j]*g;// temptation to misreport
    chi[j] = -G[j]*w;// 

//    base[j] = 1;// exp( sd_base*err_log_base[j] ); //mu_base + 
    b[j] = 2*inv_logit(rho*tau[j]) - 1;
//    R[j] = ( log1p_exp(rho*chi[j]) - log1p_exp(rho*(chi[j]-1)) )/rho;
    R[j] = inv_logit( rho*chi[j] );//+ rep_vector(sd_base*err_log_base[j], N) );
    m[j] = //base[j] * 
            b[j] .* R[j];  // extent of misreporting resulting from bias effort 
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
//    D[j,(L+1):N] = m[j,(L+1):N] - d[1]*m[j,L:(N-1)] - d[2]*m[j,(L-1):(N-2)] - d[3]*m[j,(L-2):(N-3)];     
// if estimated d < 1, this'd suggest taking more than one-period to reverse
    D[j,(L+1):N] = m[j,(L+1):N];
    for (l in 1:L) {
//      D[j,(L+1):N] += (-1)* d[l] * m[j,(L+1-l):(N-l)];
      D[j,(L+1):N] = D[j,(L+1):N] + (-1)* d[l] * m[j,(L+1-l):(N-l)];
      }
//===============================================================================

      
//===============================================================================
// Define real EM component
phi[j] = //rep_vector(p3, N) + 
          Z[j]*p;//[2:I];// temptation to manage current-period real earnings upward
    Real[j] = inv_logit(phi[j]);// ( log1p_exp(rho*phi[j]) - log1p_exp(rho*(phi[j]-1)) )/rho;
//===============================================================================


//===============================================================================
// Define the fraction of real EM that is sales-based and due to shfiting next period's sales forward
      //zeta[j] = T[j]*s;// temptation to manage current-period real earnings upward
      //sigma[j] = inv_logit(zeta[j]);//( log1p_exp(rho*zeta[j]) - log1p_exp(rho*(zeta[j]-1)) )/rho;
// fraction of Real[j] that is sales-based, rather than RnD-based,
//   where T indicates whether j is a retailer and has RnD spending or not reported in last period
//===============================================================================


//===============================================================================
// Define seasaonality component
    season_q[j] = append_row( season_raw[j], -sum(season_raw[j]) );
    season_n[j,1:N] = season_q[j,q[1:N]];
//===============================================================================


//===============================================================================
// Define real EM's LT impact on alpha
  alpha[j,1] = mu_alpha;  
  for (n in 2:N) {
    alpha[j,n] = alpha[j,n-1] - theta*Real[j,n-1]; 
//    alpha[j,n] = alpha[j,n-1] - (1 - sigma[j,n-1]) .* (theta*Real[j,n-1]); 
    }
//===============================================================================



/* ===============================================================================
    u[j,1] ~ normal(//season_n[j,1] + 
                    mu_u1, sd_y);
    u[j,2:N] ~ normal(//season_n[j,2:N] +
                          mu_alpha //alpha[j,2:N] 
                          + beta*u[j,1:(N-1)] 
#                        - sigma[j,1:(N-1)] .* Real[j,1:(N-1)]
                        , sd_y); 
===============================================================================
*/
  err_y[j,1] = r[j,1] - (season_n[j,1] + mu_u1)
                  - Real[j,1]
                  - D[j,1];

  err_y[j,2:N] = r[j,2:N] 
                  - ( season_n[j,2:N] + alpha[j,2:N] + beta*r[j,1:(N-1)] )
//                                                     - sigma[j,1:(N-1)] .* Real[j,1:(N-1)] ) 
                                             // fraction of Real[j] that is sales-based, rather than RnD-based,
                  - ( Real[j,2:N] - beta*Real[j,1:(N-1)] )
                  - ( D[j,2:N] - beta*D[j,1:(N-1)] );

    }


  // priors 

  sd_y ~ exponential(5);//2);//normal(0, 1); //sd_y ~ exponential(1); //sd_y ~ student_t(3, 0, 1); 
sd_season ~ normal(0, 0.2);//exponential(5);////0.5);//1);//exponential(1); //
mu_season ~ normal(0, 0.2);//0.5);//1);

  mu_u1 ~ normal(0, 0.5);//0.5);// //student_t(3, 0, 1);//
theta ~ normal(0.5, 0.5);//~ normal(0, 1);   exponential(2);//
//===============================================================================
//  mu_alpha ~ normal(0.5, 0.5);//normal(0, 1);//student_t(3, 0, 1);//exponential(2);//normal(0.5, 0.5);//
//  beta ~ normal(0.5, 0.5);
  ab_mu ~ normal(0.5, 0.5);//0.5);//0.25);//0.2);//normal(0, 1);//
//  ab_mu[1] ~ normal(0.5, 0.5);//0.5);//0.25);//0.2);//normal(0, 1);//
//  ab_mu[2] ~ normal(0.5, 0.5);//0.5);//0.25);//0.2);//normal(0, 1);//
  
  ab_sd ~ normal(0, 0.1);//25);//0.5);//1);//exponential(1);//
  ab_L ~ lkj_corr_cholesky(2);//lkj_corr_cholesky(2);
  ab_err ~ normal(0, 0.05);//0.1);//25);//0.5);//1); // implies:  w_raw ~ multi_normal(w_mu, quad_form_diag(w_L * w_L', w_sd));


p_raw ~ normal(0.5, 0.5); 
g_raw ~ normal(0.5, 0.5); 
w_raw ~ normal(0.5, 0.5); 
//!!!
//===============================================================================
pgw_mu ~ normal(0, 0.63);//note: half-normal(0.63) has mean around 0.5;//normal(0.5, 0.2);//exponential(2); //
  pgw_sd ~ normal(0, 0.1);//0.05);//0.2);
  pgw_L ~ lkj_corr_cholesky(2);
  pgw_err ~ normal(0, 0.05); //0.1); implies:  w_raw ~ multi_normal(w_mu, quad_form_diag(w_L * w_L', w_sd));
//===============================================================================
/*
//  g ~ normal(0, 1);//student_t(4, 0, 1);//student_t(3, 0, 5);//normal(0, 1);//10);//0); //0.2); // g_init
        //g0 ~ normal(0, 1);//2);
  g_mu[1:2] ~ normal(0, 0.5);//1);
  g_mu[3] ~ normal(0, 0.1);//1);
  g_sd ~ normal(0, 0.1);//0.2);
  g_L ~ lkj_corr_cholesky(2);
  g_err ~ normal(0, 0.1);//0.05); //0.1); implies:  w_raw ~ multi_normal(w_mu, quad_form_diag(w_L * w_L', w_sd));
//===============================================================================
//    w ~ normal(0, 0.5);//1);//student_t(4, 0, 1);//student_t(3, 0, 5);//normal(0, 5);//10);//0); //0.2); //1);  // w_init
  w_mu[1:2] ~ normal(0, 0.5);//1); //student_t(3,0,1);//4, 0, 1);//
  w_mu[3] ~ normal(0, 0.5);//1); //student_t(3,0,1);//4, 0, 1);//
  w_sd ~ normal(0, 0.1);//0.2); //2);
  w_L ~ lkj_corr_cholesky(2);
  w_err ~ normal(0, 0.1);//0.2); // implies:  w_raw ~ multi_normal(w_mu, quad_form_diag(w_L * w_L', w_sd));
//!!!    w0 ~ normal(0, 0.5);//1);//2);
//===============================================================================
*/


//  d ~ dirichlet(rep_vector(0.1, L));   // = 1.0 means the dist is uniform over the possible simplices;<1.0 toward corners 
//!!!  mu_base ~ normal(0, 1);
//  sd_base ~ exponential(3);//normal(0, 1);


  for (j in 1:J) { 

//    err_log_base[j] ~ normal(0, 1);
    season_raw[j] ~ normal(mu_season, sd_season);

    err_y[j] ~ normal(0, sd_y);



/*//!!!
//===============================================================================
    u[j,1] ~ normal(//season_n[j,1] + 
                    mu_u1, sd_y);
    u[j,2:N] ~ normal(//season_n[j,2:N] +
                          mu_alpha //alpha[j,2:N] 
                          + beta*u[j,1:(N-1)] 
#                        - sigma[j,1:(N-1)] .* Real[j,1:(N-1)]
                        , sd_y); 
//===============================================================================
*/

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

