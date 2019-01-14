functions{

  vector[] sem_mean(vector[] alpha, real[,,] B, real[,,] gamma, int[] g, int k, int Ng, int gamind, real[,] meanx){
    matrix[k,k] iden;
    vector[k] evlv[Ng];

    iden = diag_matrix(rep_vector(1.0, k));

    for(j in 1:Ng){
      if(gamind == 1){
        evlv[j] = inverse(iden - to_matrix(B[,,j])) * (alpha[j] + to_matrix(gamma[,,j]) * to_vector(meanx[,j]));

      } else {
        evlv[j] = inverse(iden - to_matrix(B[,,j])) * alpha[j];
      }
    }

    return evlv;
  }

  real sem_lv_lpdf(matrix x, real[,,] alpha, real[,,] B, real[,,] psi, real[,,] gamma, int gamind, real[,] meanx, int[] g, int k, int N, int Ng, int diagpsi, int fullbeta, int nlv, int[] lvind, int nlvno0){
    real ldetcomp[Ng];
    matrix[k,k] iden;
    vector[k] alpha2[Ng];
    vector[k] psivecinv[Ng];
    matrix[k,k] psimatinv[Ng];
    matrix[k,k] psimat[Ng];
    matrix[k,k] siginv[Ng];
    vector[k] xvec;
    vector[k] evlv[Ng];
    int idx[(k-nlv+nlvno0)];
    real xvectm;
    real ldetsum;
    int nov;
    int nidx;

    nov = k - nlv;
    nidx = nov + nlvno0;

    iden = diag_matrix(rep_vector(1.0, k));

    if(nlvno0 > 0){
      idx[1:nlvno0] = lvind;
    }
    if(nov > 0){
      for(j in 1:nov){
        idx[nlvno0+j] = nlv + j; //nlvno0 + j?
      }
    }

    for(j in 1:Ng){
      alpha2[j] = to_vector(alpha[,1,j]);
    }

    evlv = sem_mean(alpha2, B, gamma, g, k, Ng, gamind, meanx);

    if(diagpsi){
      for(j in 1:Ng){
        for(i in 1:nidx){
          psivecinv[j,idx[i]] = 1/psi[idx[i],idx[i],j];
        }
        psimatinv[j] = diag_matrix(psivecinv[j]);

        siginv[j,1:nidx,1:nidx] = (iden[idx,idx] - to_matrix(B[idx,idx,j])') * psimatinv[j,idx,idx] * (iden[idx,idx] - to_matrix(B[idx,idx,j]));

	if(fullbeta){
	  ldetcomp[j] = log_determinant(iden[idx,idx] - to_matrix(B[idx,idx,j]));
	  ldetcomp[j] = -2 * ldetcomp[j] + sum(log(diagonal(to_matrix(psi[idx,idx,j]))));
	} else {
          ldetcomp[j] = sum(log(diagonal(to_matrix(psi[idx,idx,j]))));
  	}
      }
    } else {
      for(j in 1:Ng){
	psimat[j] = to_matrix(psi[,,j]) + to_matrix(psi[,,j])' - diag_matrix(diagonal(to_matrix(psi[,,j])));

	ldetcomp[j] = log_determinant(psimat[j,idx,idx]);
	if(fullbeta){
	  ldetcomp[j] = ldetcomp[j] - 2 * log_determinant(iden[idx,idx] - to_matrix(B[idx,idx,j]));
	}

	psimatinv[j] = psimat[j];
	psimatinv[j,1:nidx,1:nidx] = inverse_spd(psimat[j,idx,idx]);
        siginv[j,1:nidx,1:nidx] = (iden[idx,idx] - to_matrix(B[idx,idx,j])') * psimatinv[j,1:nidx,1:nidx] * (iden[idx,idx] - to_matrix(B[idx,idx,j]));
      }
    }

    xvectm = 0;
    ldetsum = 0;
    for(i in 1:N){
      xvec = x[i,]';
      xvectm = xvectm + (xvec[idx] - evlv[g[i],idx])' * siginv[g[i],1:nidx,1:nidx] * (xvec[idx] - evlv[g[i],idx]);
      ldetsum = ldetsum + ldetcomp[g[i]];
    }

    return -0.5 * (ldetsum + xvectm);
  }

  matrix fill_lower(matrix x){
    matrix[rows(x),cols(x)] newx;

    newx = x;
    for(i in 1:(rows(x) - 1)){
      for(j in (i+1):rows(x)){
        newx[j,i] = x[i,j];
      }
    }
    return newx;
  }
}

data{
  int N;
  int g[N];
  int lvind[5];
  int etaind[2];
  real sampmean[5,1];
  real meanx[5,1];
  int dummyov[3];
  int dummylv[3];
  vector[2] y[N];
  vector[3] x[N];
  real lambdaframe[5,5,1];
  real thetaframe[5,5,1];
  real psiframe[5,5,1];
  real betaframe[5,5,1];
  real nuframe[5,1,1];
  real alphaframe[5,1,1];
  real lvrhoframe[5,5,1];
}

parameters{
  vector[1] lambdafree;
  vector<lower=0>[2] thetafree;
  vector<lower=0>[4] psifree;
  vector[4] betafree;
  vector[2] nufree;
  vector[3] alphafree;
  vector<lower=0,upper=1>[2] lvrhofree;
  matrix[N, 2] etavec;
}

transformed parameters{
  real lambda[5,5,1];
  real theta[5,5,1];
  matrix[2,2] thetld[1];
  real psi[5,5,1];
  real beta[5,5,1];
  real nu[5,1,1];
  real alpha[5,1,1];
  real lvrho[5,5,1];
  real mu[N,5];
  matrix[N,5] eta;

  eta = rep_matrix(0, N, 5);

  lambda = lambdaframe;
  theta = thetaframe;
  psi = psiframe;
  beta = betaframe;
  nu = nuframe;
  alpha = alphaframe;
  lvrho = lvrhoframe;

  lambda[1,1,1] = 1;
  lambda[2,1,1] = lambdafree[1];
  beta[3,2,1] = 1;
  beta[4,2,1] = betafree[1];
  beta[4,1,1] = betafree[2];
  beta[5,2,1] = betafree[3];
  beta[3,2,1] = betafree[4];
  psi[3,3,1] = 1;
  theta[1,1,1] = pow(thetafree[1],-1);
  theta[2,2,1] = pow(thetafree[2],-1);
  psi[4,4,1] = pow(psifree[1],-1);
  psi[5,5,1] = pow(psifree[2],-1);
  psi[1,1,1] = pow(psifree[3],-1);
  psi[2,2,1] = pow(psifree[4],-1);
  nu[1,1,1] = nufree[1];
  nu[2,1,1] = nufree[2];
  alpha[3,1,1] = alphafree[1];
  alpha[4,1,1] = alphafree[2];
  alpha[5,1,1] = alphafree[3];
  alpha[1,1,1] = 0;
  alpha[2,1,1] = 0;
  lvrho[1,2,1] = -1 + 2*lvrhofree[1];
  lvrho[4,5,1] = -1 + 2*lvrhofree[2];
  psi[1,2,1] = lvrho[1,2,1] * sqrt(psi[1,1,1] * psi[2,2,1]);
  psi[4,5,1] = lvrho[4,5,1] * sqrt(psi[4,4,1] * psi[5,5,1]);

  // mu definitions
  for(i in 1:N) {
    eta[i,1:2] = etavec[i];
    eta[i,3:5] = x[i]';

    mu[i,1] = nu[1,1,g[i]] + lambda[1,1,g[i]]*eta[i,1];
    mu[i,2] = nu[2,1,g[i]] + lambda[2,1,g[i]]*eta[i,1];
  }

  for(j in 1:1){
    thetld[j,1,1] = theta[1,1,j];
    thetld[j,1,2] = theta[1,2,j];
    thetld[j,2,2] = theta[2,2,j];
    thetld[j] = fill_lower(thetld[j]);
    thetld[j] = cholesky_decompose(thetld[j]);
  }

}

model {
  for(i in 1:N) {
    y[i] ~ multi_normal_cholesky(to_vector(mu[i,1:2]), thetld[g[i]]);
  }

  eta ~ sem_lv_lpdf(alpha, beta, psi, beta, 0, meanx, g, 5, N, 1, 0, 0, 2, etaind, 2);

  // Priors
  lambdafree[1] ~ normal(0,10);
  betafree[1] ~ normal(0,10);
  betafree[2] ~ normal(0,10);
  betafree[3] ~ normal(0,10);
  betafree[4] ~ normal(0,10);
  thetafree[1] ~ gamma(1,.5);
  thetafree[2] ~ gamma(1,.5);
  psifree[1] ~ gamma(1,.5);
  psifree[2] ~ gamma(1,.5);
  psifree[3] ~ gamma(1,.5);
  psifree[4] ~ gamma(1,.5);
  nufree[1] ~ normal(0,1000^.5);
  nufree[2] ~ normal(0,1000^.5);
  alphafree[1] ~ normal(0,1000^.5);
  alphafree[2] ~ normal(0,1000^.5);
  alphafree[3] ~ normal(0,1000^.5);
  lvrhofree[1] ~ beta(1,1);
  lvrhofree[2] ~ beta(1,1);
}