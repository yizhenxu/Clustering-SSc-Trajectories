scode = "
data {
  int<lower=0> Nsubs; //number of subjects

  int<lower=0> NpredsX; // ncol of design matrix X
  
  int<lower=0> NpredsS; // ncol of design matrix S

  int<lower=0> Nobs; //number of non-NA observations
  int<lower=1, upper=Nsubs> subject[Nobs]; // created fake ID 1,...,Nsubs

  matrix[Nobs,2] Y; //  outcomes matrix 
  matrix[Nobs,NpredsX] X; // design matrix X cluster-specific
  matrix[Nobs,NpredsS] S; // design matrix S shared by clusters

  int<lower=0> L; // number of mixture components, 2
  int<lower=0> nb; //number of random effects, 2

  real Gnu; // hyperparameters for sampling G
  matrix[2,2] GS;

  real Rnu; // hyperparameters for sampling R
  matrix[2,2] RS;

  int lsub[Nsubs]; //element i being total number of time points for person i
  int endidx[Nsubs]; // cumulative_sum(lsub);

  matrix[Nsubs, 5] Z;
}
parameters {
  cov_matrix[2] R; // variance of error term epsilon

  matrix[nb, Nsubs] z_1;
  matrix[2, NpredsX] beta_1;
  cov_matrix[nb] G_1; // covariance of b's

  matrix[nb, Nsubs] z_2;
  matrix[2, NpredsX] beta_2;
  cov_matrix[nb] G_2;

  vector[NpredsS] betaS[2]; //beta for S
  
  vector[5] alpha; 
}
transformed parameters {
  
  simplex[2] theta[Nsubs];
  
  matrix[Nsubs, nb] b_1; // b0(Y1), b1(Y1), b0(Y2), b1(Y2)
  matrix[nb,nb] sqrG_1;
  matrix[Nobs,2] mu_1; // X * beta + Z * b

  matrix[Nsubs, nb] b_2;
  matrix[nb,nb] sqrG_2;
  matrix[Nobs,2] mu_2;

  sqrG_1 = cholesky_decompose(G_1);
  b_1 = (sqrG_1 * z_1)'; // b[n] ~ N(0, L*L') = N(0, G)
  for(n in 1:Nobs){
    mu_1[n,1] = S[n] * betaS[1] + X[n] * beta_1[1]' + b_1[subject[n],1];
    mu_1[n,2] = S[n] * betaS[2] + X[n] * beta_1[2]' + b_1[subject[n],2];
  }

  sqrG_2 = cholesky_decompose(G_2);
  b_2 = (sqrG_2 * z_2)';
  for(n in 1:Nobs){
    mu_2[n,1] =  S[n] * betaS[1] + X[n] * beta_2[1]' + b_2[subject[n],1];
    mu_2[n,2] =  S[n] * betaS[2] + X[n] * beta_2[2]' + b_2[subject[n],2];
  }
  
  for(i in 1:Nsubs){
    theta[i,1] = inv_logit(Z[i] * alpha);
    theta[i,2] = 1 - theta[i,1];
  }

}
model {

  to_vector(z_1) ~ normal(0,1);
  to_vector(z_2) ~ normal(0,1);

  betaS[1] ~ normal(0,10);
  betaS[2] ~ normal(0,10);
  alpha ~ normal(0,10);
  
  to_vector(beta_1) ~ normal(0,10); // Gaussian(0, var = 100)
  G_1 ~ inv_wishart(Gnu, GS);

  to_vector(beta_2) ~ normal(0,10);
  G_2 ~ inv_wishart(Gnu, GS);

  R ~ inv_wishart(Rnu, RS);


  for(i in 1:Nsubs){
    vector[L] lps = log(theta[i]);
      for(n in (endidx[i]-lsub[i]+1):endidx[i] ){
        lps[1] += multi_normal_lpdf(Y[n] | mu_1[n], R);
        lps[2] += multi_normal_lpdf(Y[n] | mu_2[n], R);
      }//Ni
      target += log_sum_exp(lps);
  }//i

}
generated quantities {
  matrix[Nsubs, L] logpyz;
  real log_lik = 0;

  for(i in 1:Nsubs){
    vector[L] lps = log(theta[i]);
      for(n in (endidx[i]-lsub[i]+1):endidx[i] ){
        lps[1] += multi_normal_lpdf(Y[n] | mu_1[n], R);
        lps[2] += multi_normal_lpdf(Y[n] | mu_2[n], R);
      }//Ni
      for(l in 1:L){
        logpyz[i,l] = lps[l];
      }//L
    log_lik += log_sum_exp(lps);
  }//i


}
"


datalist = list(Nsubs = max(refdat$ID),
                Nobs = nrow(xY), NpredsX = ncol(xY), NpredsS = ncol(xS),
                subject = refdat$ID, Y = Yvec,
                X = xY, S = xS,
                nb = 2, L = 2, Gnu = 4, Rnu = 4, GS = diag(2), RS = diag(2),
                lsub = lsub, endidx = cumsum(lsub),
                Z = Z)


parvec = c("betaS","beta_1","beta_2","G_1","G_2","b_1","b_2","R", "alpha", "log_lik", "logpyz")


