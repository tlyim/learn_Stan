data {
  int<lower=1> N;
  real x[N];

  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
}

transformed data {
/* Note: added a nugget or jitter of 10^−10 so the marginal covariance matrix before computing its Cholesky decomposition in order to stabilize the numerical calculations. Often the nugget is taking to be square root of the floating point precision, which would be 10−8 for double precision calculations.
 https://betanalpha.github.io/assets/case_studies/gp_part1/part1.html
*/
  matrix[N, N] cov = cov_exp_quad(x, alpha, rho) + diag_matrix(rep_vector(1e-10, N));
  matrix[N, N] L_cov = cholesky_decompose(cov);
}

parameters {}
model {}

generated quantities {
  vector[N] f = multi_normal_cholesky_rng(rep_vector(0, N), L_cov);
  vector[N] y;
  for (n in 1:N)
    y[n] = normal_rng(f[n], sigma);
}



