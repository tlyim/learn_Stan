data {
  int<lower=0> J; // number of firms
  int<lower=0> N; // number of observations (for different years)
  int<lower=0> Q; // number of seasons, eg, 4 for quarters
  int<lower=1,upper=Q> q[N]; // index for seasons, eg, q[n] = 1 if n belongs to Q1
  int<lower=0> I; // number of coefficents for the column vectors in covariate matrix Z (which can equal X)
  int<lower=0> K; // number of coefficents for the column vectors in covariate matrix X
  int<lower=0> H; // number of coefficents for the column vectors in covariate matrix G
  int<lower=0> L; // length of simplex vector representing reversal lasting L periods

  vector[N] r[J]; // reported figures (scaled by Rev, total revenue)
  matrix[N,I] Z[J]; // covariates (capturing the dis/incentives of real earnings management)
  matrix[N,K] X[J]; // covariates (capturing the dis/incentives of misreporting)
  matrix[N,H] G[J]; // covariates (capturing the strength of internal governance and external monitoring mechansims)
}
transformed data {
}
parameters {
  real mu_u1;
  real mu_alpha; // intercept coefficient (drift) of the AR(1) process of the unbiased figure y[n]
  real<lower=0,upper=1> beta; // slope coefficient of the AR(1) process of y[n]
  real<lower=0> sd_y; // sd of the underlying unbiased figures (vector y)

}
transformed parameters {
//vector[N] u[J];
vector[N] err[J];  // error term in p
  
  for (j in 1:J) {  
    err[j,1] = r[j,1]- mu_u1; //season_n[j,1] + 
//=======================================================      
    for (n in 2:N) { 
      err[j,n] = r[j,n] - mu_alpha + beta*r[j,n-1];//season_n[j,n] + 
//      u[j,n] = alpha[j,n-1] + beta*u[j,n-1] + sd_y*err[j,n];//season_n[j,n] + 
//      u[j,n] = season_n[j,n] + mu_alpha + beta*u[j,n-1];
    }

//    u[j] = r[j];// - D[j];// - RealST[j];
    
  }
  
}
model {

// priors 
mu_u1 ~ normal(0, 1); //student_t(3, 0, 1);//
mu_alpha ~ normal(0, 1);//student_t(3, 0, 1);//

// priors for non-negative, eg, sd
beta ~ normal(0.5, 0.5);
sd_y ~ exponential(1); 
//sd_y ~ student_t(3, 0, 1); 
//sd_y ~ normal(0, 1); 


  for (j in 1:J) { 
    err[j] ~ normal(0, sd_y);
//    r[j,1] ~ normal(mu_u1, sd_y);
//    for (n in 2:N) { 
//      r[j,n] ~ normal(mu_alpha + beta*r[j,n-1], sd_y);//season_n[j,n] + 
//      }    
    }

}
generated quantities {
}

