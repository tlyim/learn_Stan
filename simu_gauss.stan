data {
  int<lower=1> N;
  real x[N];

  real<lower=0> alpha;
  real<lower=0> rho;
  real<lower=0> sigma;
}


// transformed data{} means the way to do simulation (?)
// rep_vector(1e-10, N) to create a n-sized vector of 1e-10 (?), with diag_matrix() to turn it into a diagonal matrix (?) 
// The Cholesky decomposition of a real-valued symmetric positive-definite matrix A is a decomposition of the form: A = LL^T, 
// where L is a lower triangular matrix with real and positive diagonal entries, and L^T denotes the conjugate transpose of L. L is unique if A is positive-definite, not just semi-positive-definite.
// The Cholesky decomposition is roughly twice as efficient as the LU decomposition for solving systems of linear equations.
// add a nugget or jitter of 10^−10 to the marginal covariance matrix before computing its Cholesky decomposition in order to stabilize the numerical calculations
// Often the nugget is taking to be square root of the floating point precision, which would be 10−8 for double precision calculations.
transformed data {
  matrix[N, N] cov =   cov_exp_quad(x, alpha, rho)
                     + diag_matrix(rep_vector(1e-10, N));
  matrix[N, N] L_cov = cholesky_decompose(cov);
}
// (Note: cov is the A and L_cov is the L of the Cholesky decomposition of A)


// Simulation, not estimation; that's why null parameters{} and model{} to fit (?)
parameters {}
model {}

// _rng() in generated quantities{} means the step to actually simulate the data (?)
// rep_vector(0, N) to create a n-sized vector of 0 (?), which is the mean function of the simulated GP (?)
// multi_normal_cholesky_ indicates that the multi-normal distribution simulation of f(x) is based on the L_cov of the Cholesky decomposition, rather than based on the original covariance matrix cov
// for (n in 1:N)_; loops over to generate the N observations of y based of the N simulated f(x) (?)
generated quantities {
  vector[N] f = multi_normal_cholesky_rng(rep_vector(0, N), L_cov);
  vector[N] y;
  for (n in 1:N)
    y[n] = normal_rng(f[n], sigma);
}
