#' ---
#' title: "State space model of accounting manipulation - Stan"
#' output: html_notebook
#' fig_width: 6 
#' fig_height: 6 
#' ---
#' 
## ------------------------------------------------------------------------

#knitr::purl("ssm_Stan0.5.1dev.Rmd", output = "ssm_Stan0.5.1dev.R", documentation = 2)
.libPaths()
.libPaths( c( "~/R_library/3.5", .libPaths() ) )
#.libPaths( c( .libPaths(), "~/R_library/3.5" ) )
.libPaths()

#library(tidyverse)
library(dplyr)
library(magrittr)

#load previously saved Stan model objects
#load(file = "StanObjects.rda")

#library(shinystan)
library(bayesplot)
library(rstan)
rstan_options(auto_write = TRUE) # avoid recompilation of unchanged Stan programs

# The following can affect distributional computing over a cluster
options(mc.cores = parallel::detectCores())  # 

# The following throws an error in compiling .stan
#Sys.setenv(LOCAL_CPPFLAGS = '-march=native') # For improved execution time but can throw errors for some CPUs


#' 
#' 
#' # Estimate the simulated model with Stan - misreporting extent 
#'   (forecasts specified in generated quantities block)
#' //ignore these lines: 
#' //  COGS2Rev, SGA2Rev, RnD2Rev with Rev2Rev = unit vector as the intercept)
#' //   int<lower=0,upper=1> M[N];   // M = 1 for the decision to misreport; = 0 if report honestly  
#' 
## NA
#' 
#' # Simulate data
## NA
#' 
#' 
#' 
#' # Simulate data and fit the debug model (also, save Stan model objects) 
## ------------------------------------------------------------------------
# J=80,N=4*30: Apparently, no need to go beyond this (take WRDS about 2.5 hrs to finish)
# J=45,N=4*30 gives reasonable estimates of the baseline model of AR(1) with hierarchical seasonability
# J=60,N=4*20 also seems to be acceptable
# Can consider J=42,N=4*18 for model dev process only  
#   but J=45,N=4*20 gives more reasonable estimates  
#===========================================
J = 80#42#40#45#20 #200     # number of firms
N = 4*25#30#18#15#5#15   #number of (quarterly) observations of a firm
Q = 4 # quarterly time series

rho = 1#2#5  # curvature of the soft-clipping function; rho=50 would make nearly hard-clipping
mu_base = -1.6
sd_base = 0.2 #0.05  
mu_alpha = 0.04 #0.1#0.2  # intercept parameter of the AR(1) processs of the underlying unbiased figure
#mu_alpha = 0.04  # intercept parameter of the AR(1) processs of the underlying unbiased figure
beta = 0.6#0.85  # slope parameter of the AR(1) processs of the underlying unbiased figure
theta = 0.15 #0.2 #0.1 # real EM's LT impact on alpha[j,n]
y_LT = mu_alpha/(1-beta)
mu_u1 = y_LT # initial y for data simulation set to three times of the LT stationary level
ab_mu = c(mu_alpha, beta)
sd_y = 0.08  # sd of underlying unbiased figures is likely to be industry-specific / business model specific
mu_season = c(-0.12, -0.06, 0.15)
sd_season = 0.1

  
  #delta = 0.65 # reversal coeff. of real EM (eg, prior year's heavily discounted sales taking from this year's)
  # delta can be bigger than, eg, sales taken from last year's heavy discount that would have been higher w/o discount; but delta can be less than 1 if the heavy discount does stimulate more sales than taken from next yaer's

#pvec <- matrix(data = p, nrow = I)  # define column vector as a matrix[N,1]
#Z.l[[3]]%*%p %>%  str()  #ok with p as will be converted accordingly


# parameters of the DGP of the likelihood of using sales- vs. RnD-based real EM (sigma = 1 or 0)
s = c(0, 0.75, 0.18)# # soft-clipping  rho=5
#s = c(0.4, 0.7, 0.5)# # soft-clipping  rho=5
S = length(s)

# parameters of the data generation process of the temptation to do real EM
#p = c(0.95, 0.7)# middle: 3.0, 
#!!! p = c(-0.22, 3.0, 0.13) # exp(p[1]) = 0.802
#p = c(3.0, 0.13) 
p = c(0.62, 0.45, 0.3) #
#p = c(0.7, 0.3, 0.13) #[before rescaling Z2] c(0.7, 3.0, 0.13)
I = length(p) # soft-clipping  rho=5

# parameters of the data generation process of the temptation to misreport
#g = c(0.77, 0.5) #c(0.5, 0.5)
g = c(0.3, 0.2, 0.4) 
#g = c(0, 0.11, 0.22) #[before rescaled X2] g = c(0.03, 2.3)   # g = c(0.5, 1)
K = length(g)

# parameters of the data generation process of the strength of governance and monitoring 
# (eg, audit committee chair with financial expertise, big 4 auditor, etc)
#w = c(0.1, 0.57) #w = c(0, 0)
    #w = c(0, 0.43, 0.85) 
w = c(0.72, 0.15, 0.62)  
#[old] w = c(3.7, 4.3, 0.6)
H = length(w)

d = c(0.05, 0.7, 0.25) # fractional reversal of prior-period manipluation by accruals
L = length(d)

NoReportedRnD_freq = 0.85
IsRetailer_freq = 0.7
ratioNED_req = 0.25 #0.3  # min ratio of non-executive directors (NED) required by law
Big4_freq = 0.7#0.6
downgrade_freq = 0.3 #0.5
integrity = 0.4 #0.35 #(0,1)  0.5 = average       
# integrity level, between (0,1), of the society, interpreted as likelihood to behave opportunistically
# integrity level of an individual, indiv_integ between (-Inf,Inf), is 
#   affected by a society's general integrity level in a stochastic fashion varying from period to period

#M <- rbinom(N, 1, integrity)   # the decision to misreport earnings

# Define function to convert data.frame with gvkey as the key to a list of matices
df2mat.l <- function(z){
  z %>% 
  split(.$gvkey) %>% 
  lapply(function(z) select(z, -c(gvkey))) %>% 
  lapply(as.matrix)
  }

#set.seed(8)
# T.l captures firm characteristics related to whether firms are more likely to use sales- vs. RnD-based real EM
T.l <- data.frame(
  gvkey = rep(1:J, times=N),                   
  constant = rep(1, J*N),
  IsRetailer = rbinom(J*N, 1, IsRetailer_freq),
  NoReportedRnD = (10)*rbinom(J*N, 1, NoReportedRnD_freq)#,
 ) %>%
  df2mat.l()

# make Z.l different from X.l to avoid less than full rank in estimation
Z.l <- data.frame(
  gvkey = rep(1:J, times=N),                   
  constant = (-10)*rep(1, J*N),
  distortion = (10)*rbeta(J*N, 3*(1-integrity), 3*integrity),#*rbeta(J*N, integrity, 1-integrity), 
  downgrade = rbinom(J*N, 1, downgrade_freq)#,
#constant = rep(1, J*N)#,
 ) %>%
  df2mat.l()

# covariate matrix X capturing the dis/incentive to misreport
# downgrade = 1 if one or more analysts following the firm downgraded their stock recommendation in the last period
X.l <- data.frame(
  gvkey = rep(1:J, times=N),                   
  constant = rep(1, J*N),
  aggressive = (10)*( rbeta(J*N, 3*(1-integrity), 3*integrity) - 0.5 ), 
  downgrade = (10)*rbinom(J*N, 1, downgrade_freq)#,
  ) %>%
  df2mat.l()

get_ratioNED <- function(ratioNED_req, J, N) {
  CC = pmax(ratioNED_req, rbinom(J, 20, ratioNED_req)/20) 
  EE = rep(NULL, J)
  for (n in 1:N)  { 
    A = rbinom(J, 1, 0.5)/20
    B = -rbinom(J, 1, 0.5)/20
    CC = pmin(0.65, pmax(ratioNED_req, A+B+CC)) #%>%  print()
    EE = c(EE, CC)
    } 
  return(EE)
  }
# covariate matrix G capturing the strength of internal governance and external monitoring mechansims
G.l <- data.frame(
  gvkey = rep(1:J, times=N), # note the difference from times=N
  constant = (1)*rep(1, J*N),
  ratioNED = (10)*get_ratioNED(ratioNED_req, J, N), 
  Big4 = (2)*rbinom(J*N, 1, Big4_freq)#,
  ) %>%
  df2mat.l()
# G.l %>%  arrange(gvkey) %>% head(120)
#colnames(Z) <- c("Rev2Rev", "COGS2Rev", "SGA2Rev", "RnD2Rev")

simu_SSM_data <- list(J = J, 
                      N = N, 
                      Q = Q,
                      S = S,
                      I = I,
                      K = K,
                      H = H,
                      L = L,
                      T = T.l,
                      Z = Z.l,
                      X = X.l,
                      G = G.l,
rho = rho,
                      sd_y = sd_y,
                      mu_season = mu_season,
                      sd_season = sd_season,
                      mu_u1 = mu_u1,
                      theta = theta,
                      mu_alpha = mu_alpha,
                      beta = beta, 
#sd_alpha = sd_alpha,
#                      sd_gamma = sd_gamma,
#                      sd_omega = sd_omega,
                      mu_base = mu_base,
                      sd_base = sd_base,
                      s = s,
                      p = p,
                      g = g,
                      w = w,
                      d = d
                      )
# save(simu_SSM_data, file="simu_SSM_data.rda")
# load(file="simu_SSM_data.rda")

# Simulate SSM data
simu_fit <- stan(#model_code = simu_SSM@model_code,
                 file = "simu_SSM0.5.1dev.stan",  # Stan program
#                 file = "simu_SSM2.4.stan",  # Stan program
                 data=simu_SSM_data, iter=1, #seed=987,  
                 chains=1, algorithm="Fixed_param")  

q <- extract(simu_fit)$q[1,] 

r.mat <- extract(simu_fit)$r[1,,]
y.mat <- extract(simu_fit)$y[1,,]
Real.mat <- extract(simu_fit)$Real[1,,]
m.mat <- extract(simu_fit)$m[1,,]
D.mat <- extract(simu_fit)$D[1,,]
b.mat <- extract(simu_fit)$b[1,,]
R.mat <- extract(simu_fit)$R[1,,]
sigma.mat <- extract(simu_fit)$sigma[1,,]

Real.mat[1,] %>% round(4)
Real.mat %>% round(4) %>%  mean()
Real.mat %>% round(4) %>%  median()
Real.mat %>% round(4) %>%  max()
Real.mat %>% round(4) %>%  min()
cat("\n\n")
r.mat[1,] %>% round(4)
y.mat[1,] %>% round(4)
y.mat[,1] %>% round(4) %>%  mean()
y.mat[,1] %>% round(4) %>%  max()
y.mat[,1] %>% round(4) %>%  min()

b.mat[1,] %>% round(4)
b.mat %>% round(4) %>%  max()
b.mat %>% round(4) %>%  min()

R.mat[1,] %>% round(4)
R.mat %>% round(4) %>%  max()
R.mat %>% round(4) %>%  min()
cat("\n\n")
m.mat[1,] %>% round(4)
D.mat[1,] %>% round(4)
D.mat %>% round(4) %>%  mean()
D.mat %>% round(4) %>%  max()
D.mat %>% round(4) %>%  min()
cat("\n\n")
sigma.mat[1,] %>% round(4)

# Prepare SSM data for estimation using MCMC
#N_new = 3
SSM_data <- list(J = J,
                 N = N, 
                 Q = Q,
                 q = q,   
                 S = S,
                 I = I,
                 K = K,
                 H = H,
                 L = L,
#mu_base = mu_base,
rho = rho,
                 r = r.mat, 
                 T = T.l,
                 Z = Z.l,
                 X = X.l,
                 G = G.l#,
                 # Z_new = Z[1:N_new,],
                 # N_new = N_new
                 )


#' 
#' # Fit the full model (by referring to the debug model) to estimate underlying parameters
## ------------------------------------------------------------------------

n_iter = 2000#3000
n_refresh = max(2, n_iter/10)
n_warmup = 1000
#n_warmup = n_iter/2 
n_chains = 6



initf2 <- function(chain_id = 1) {
  list(mu_alpha = mu_alpha, beta = beta, theta = theta, ab_mu = ab_mu, #shape = shape,
       s = s, p = p, g = g, w = w, d = d, mu_season = mu_season, sd_season = sd_season,
       mu_u1 = mu_u1, #u = array(y_LT, dim = c(J,N)),
       rho = rho, sd_y = sd_y, mu_base = mu_base, sd_base = sd_base) #, sd_gamma = sd_gamma, sd_omega = sd_omega)
  }
# generate a list of lists to specify initial values
init_ll <- lapply(1:n_chains, function(id) initf2(chain_id = id))

# Run the full model and refer to the debug model to save compilation time 
system.time({
fit2 <- stan(
#fit1 <- stan(
  
file = "SSM0.5.1dev.stan",  # Stan program
#  file = "SSM2.4.stan",  # Stan program
  data = SSM_data,    # named list of data
#  model_code = SSM@model_code, 
#  fit = fit0_debug,   # to save compilation time if the debug model was run
  control = list(adapt_delta = 0.8,#9, 
#                 stepsize = 0.05, #0.05, #0.01
                 max_treedepth = 15),#15),    # adjust when there're divergent transitions after warmup
  init = init_ll, 
  chains = n_chains,#1,#8,             # default is 4 Markov chains
#  cores = 8,
#  seed = 555, 
  iter = n_iter,#4000,#7000,#12000, #500 #2000,            # total number of iterations per chain
  warmup = n_warmup,#500,          # number of warmup iterations per chain
  refresh = n_refresh#500          # = 0 means no progress shown
  )
})

#fit2 <- fit1
#save(fit2 , file="fit2.rda")#, file="/scratch/city/fit2.rda")
# beepr::beep()

#extract(fit1)$y[1,]
#pairs(extract(fit1), pars = c(a, b)) #, y[1], sd_y, m[1], sd_m, c[1], c[2], c[3]))


#' 
#' # Summary of the posterior sample 
## ------------------------------------------------------------------------
#====================================================
  # sd_y = 0.08, mu_u1 = y_LT = 0.1; mu_alpha = 0.04; beta = 0.6;  
# theta = 0.15; rho = 1, # sd_base = 0.2; mu_base = -1.6 
  # sd_season = 0.1, mu_season = c(-0.12, -0.06, 0.15); 
# s = c(0, 0.75, 0.18)   
  # p = c(0.62, 0.45, 0.3) 
  # g = c(0.3, 0.2, 0.4) 
  # w = c(0.72, 0.15, 0.62) 
# d = c(0.05, 0.7, 0.25) 
#====================================================
print(fit2, c(
              "sd_y", "mu_u1", "mu_alpha", "beta", #"rho", 
              "theta", "sd_season", "mu_season", #"sd_base", #"mu_base", 
              "s",
              "p",
              "g", 
              "w", 
              "d", 
              "gw_mu", "gw_sd", "gw_L", "gw_err", 
              "p_mu", "p_sd", "p_L", "p_err", 
              "ab_mu", "ab_sd", "ab_L", "ab_err"
              ), 
      probs = c(0.1, 0.5, 0.9), digits = 3) #, probs = c(.05,.95))


#' 
#' 
#' # Pairs plots
## ------------------------------------------------------------------------
posterior <- as.array(fit2)
lp <- log_posterior(fit2)
np <- nuts_params(fit2)  # nuts parameters
#par(mar = c(4, 4, 0.5, 0.5))
knitr::opts_chunk$set(fig.width=6, fig.height=6)
color_scheme_set("gray") #"darkgray")
#Error: cannot allocate vector of size 782.8 Mb
#mcmc_parcoord(posterior, np = np)

#http://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html
#The main problem is that large steps are required to explore the less narrow regions efficiently, but those steps become too large for navigating the narrow region. The required step size is connected to the value of tau. When tau is large it allows for large variation in theta (and requires large steps) while small tau requires small steps in theta.
#The non-centered parameterization avoids this by sampling the eta parameter which, unlike theta, is a priori independent of tau. Then theta is computed deterministically from the parameters eta, mu and tau afterwards. Here’s the same plot as above, but with eta[1] from non-centered parameterization instead of theta[1] from the centered parameterization:

mcmc_pairs(posterior, np = np, pars = c("mu_u1", "sd_y", "s[1]", "s[2]", "s[3]"))
mcmc_pairs(posterior, np = np, pars = c("ab_mu[1]", "ab_mu[2]", "s[1]", "s[2]", "s[3]"))
mcmc_pairs(posterior, np = np, pars = c("theta", "mu_season[1]", "s[1]", "s[2]", "s[3]"))
mcmc_pairs(posterior, np = np, pars = c("mu_season[2]", "mu_season[3]", "s[1]", "s[2]", "s[3]"))
mcmc_pairs(posterior, np = np, pars = c("gw_mu[1]", "gw_mu[2]", "s[1]", "s[2]", "s[3]"))
mcmc_pairs(posterior, np = np, pars = c("gw_mu[3]", "gw_mu[4]", "s[1]", "s[2]", "s[3]"))
mcmc_pairs(posterior, np = np, pars = c("gw_mu[5]", "gw_mu[6]", "s[1]", "s[2]", "s[3]"))
mcmc_pairs(posterior, np = np, pars = c("p_mu[1]", "p_mu[2]", "s[1]", "s[2]", "s[3]"))
mcmc_pairs(posterior, np = np, pars = c("p_mu[3]", "d[1]", "s[1]", "s[2]", "s[3]"))
mcmc_pairs(posterior, np = np, pars = c("d[2]", "d[3]", "s[1]", "s[2]", "s[3]"))

#====================================================
  # mu_u1 = y_LT = 0.1; mu_alpha = 0.04; beta = 0.6;  
# theta = 0.15; rho = 1,  
  # sd_base = 0.2; mu_base = -1.6 
  # sd_y = 0.08, sd_season = 0.1, mu_season = c(-0.12, -0.06, 0.15); 
# s = c(0, 0.75, 0.18)   
  # p = c(0.62, 0.45, 0.3) 
  # g = c(0.3, 0.2, 0.4) 
  # w = c(0.72, 0.15, 0.62) 
# d = c(0.05, 0.7, 0.25) 

