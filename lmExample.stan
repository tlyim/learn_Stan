// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; // number of observations (for different CEOs/firms)
  int<lower=0> K; // number of coefficents (/predictors)
  vector[N] y; // dependent variable
  matrix[N,K] X; // predictor variables
  int<lower=0,upper=1> M[N];   // M = 1 for the decision to misreport; = 0 if report honestly  

// forecasts  
  int<lower=0> N_new; // number of predictions
  matrix[N_new,K] X_new; // 
  
// GP data
  int<lower=0> G;
  real z[G];     // ie, x of the GP cov matrix
  vector[G] m;  // misreporting extent
}
transformed data {
// GP data
}
parameters {
  vector[K] b; // coefficients of the predictor variables
  real<lower=0> sd_y; // sd of error term
  real<lower=0,upper=1> integrity;     // integrity level of the society affecting the decision to misreport or not
//  vector[N] integrity;     // integrity level of each CEO affecting the decision to misreport or not

// GP parameters
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}
transformed parameters {
  real ilogodds;     // integrity level of the society affecting the decision to misreport or not
  ilogodds = logit(integrity);   
}
model {
  matrix[G, G] Kcov;
  matrix[G, G] Lcov;

// b ~ cauchy(0, 2.5); // common prior for each b[K]
  y ~ normal(X*b, sd_y);
  M ~ bernoulli_logit(ilogodds); 

  rho ~ gamma(2, 20);
  alpha ~ normal(0, 1);
  sigma ~ normal(0, 1);

  Kcov = cov_exp_quad(z, alpha, rho) + diag_matrix(rep_vector(square(sigma), G));
  Lcov = cholesky_decompose(Kcov);

  m ~ multi_normal_cholesky(rep_vector(0, G), Lcov);   // rep_vector(0, G) is the zero mean assumed
}

generated quantities {
  vector[N_new] y_new;
  for (n in 1:N_new)
    y_new[n] = normal_rng(X_new[n] * b, sd_y);
}


/*
   CEO's choice of misreporting (explore) vs. truth-telling (exploit) is modelled as discrete choice
     use bernoulli_logit() ?
   The stochastic process (GP? with -ve due to accrual reversal constraint?) of maximally feasible extent of misreporting
     is exogenous to the CEO but the realized value of this in a period affects the achievable reward of misreporting
     GP mean depends on (slowly changing) time-varying firm-specific, industry-specific, and macroeconomic factors ?
   Assume for simplicity that conditional on the decision to misreport, it is always optimal to choose the max extent
     Determinants of explore-exploit include 
       parameters: CEO-specific Integrity, 
       variables: firm-specific ERC, allowable misreporting limit
   Misreported CEO also learns directly the updated ERC for his firm, the incremental effect of misreporting, and
     *most importantly* the likelihood of being caught through the following period (because he knows there's misreporting). 
     Truth-telling would be harder to learn the likelihood of being caught because observing no restatement by other firms
       cannot clearly tell whether there was no misreporting or misreporting remains not caught. 
     Investors are smart enough to condition their ERC response according to the firm characteristics, etc
       For example, in bad economic times, the reward from misreporting becomes relatively higher. 
    Background noise (ie, likelihood or errors as opposed to irregularities) affects the incentive to misreport   
   To begin, assume CEO works for only one firm for his whole life. 
   Keep in mind that the unrestated earnings is a mix of misreported and truth-told earnings
      Also, assume
      - an exogenous fraction of all misreporting will be discovered in the next period 
          (regardless of the misreporting extent in an instance)
      - the fraction follows a stohastic process (again, with the mean slowly changing over time)
   Thus, the count of restatement each period, adjusted by a multiple, is an estimate of accumulated misreporting
     instances in the past     
*/

