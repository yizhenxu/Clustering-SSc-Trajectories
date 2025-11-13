
## create covariance matrix for random error terms

create_sigma2 <- function(datlist, R2){

  rowlist <- sapply(datlist, nrow)

  # no observed value for all measures
  if (all(rowlist == 0)){

    Sigma2 <- NULL

  }else{

    # if there is no observed value for any one of the measures, remove empty dataframe
    if (any(rowlist == 0)){

      datlist <- datlist[rowlist != 0]

    }

    Sigma2 <- NULL

    for (i in 1:length(datlist)){

      Cmat <- NULL

      for (j in 1:length(datlist)){

        time1 <- datlist[[i]][, 2]; time2 <- datlist[[j]][, 2]

        if (i == j){

          # diagonal matrices
          cmat <- R2[i, j] * diag(1, length(time1))

        }else{

          # off diagonal matrices
          covmat <- outer(1:length(time1), 1:length(time2), FUN = "paste", sep = ",")

          comb <- cbind(1:length(time1), match(time1, time2, nomatch = NA))


          if(all(is.na(comb[, 2]))){

            cmat <- matrix(0, nrow = length(time1), ncol = length(time2))

          }else{

            covmat[covmat %in% apply(comb, 1, paste, collapse = "," )] <- 1

            covmat[covmat != 1] <- 0

            cmat <- (covmat %>% as.data.frame %>%
                       sapply(function(x){as.numeric(levels(as.factor(x))[as.factor(x)])})) * R2[i, j]
          }
        }

        if(is.null(dim(cmat))){

          Cmat <- matrix(c(Cmat, cmat), nrow = 1)

        }else{

          Cmat <- cbind(Cmat, cmat)
        }

      }

      Sigma2 <- rbind(Sigma2, Cmat)

    }

  }

  return(Sigma2)

}


create_sigma2_b0 <- function(datlist, R2){

  rowlist <- sapply(datlist, nrow)

  # no observed value for all measures
  if (all(rowlist == 0)){

    Sigma2 <- NULL

  }else{

    # if there is no observed value for any one of the measures, remove empty dataframe
    if (any(rowlist == 0)){

      datlist <- datlist[rowlist != 0]

    }

    Sigma2 <- NULL

    for (i in 1:length(datlist)){

      Cmat <- NULL

      for (j in 1:length(datlist)){

        time1 <- datlist[[i]][, 2]; time2 <- datlist[[j]][, 2]

        if (i == j){

          # diagonal matrices
          cmat <- R2[i, j] * diag(1, length(time1))

        }else{

          # off diagonal matrices
          covmat <- outer(1:length(time1), 1:length(time2), FUN = "paste", sep = ",")

          comb <- cbind(1:length(time1), match(time1, time2, nomatch = NA))


          if(all(is.na(comb[, 2]))){

            cmat <- matrix(0, nrow = length(time1), ncol = length(time2))

          }else{

            covmat[covmat %in% apply(comb, 1, paste, collapse = "," )] <- 1

            covmat[covmat != 1] <- 0

            cmat <- (covmat %>% as.data.frame %>%
                       sapply(function(x){as.numeric(levels(as.factor(x))[as.factor(x)])})) * R2[i, j]
          }
        }

        if(is.null(dim(cmat))){

          Cmat <- matrix(c(Cmat, cmat), nrow = 1)

        }else{

          Cmat <- cbind(Cmat, cmat)
        }

      }

      Sigma2 <- rbind(Sigma2, Cmat)

    }

  }

  return(Sigma2)

}


calcpostzb0bt = function(Xi, Zi, Yi, beta1j, beta2j, G1j, G2j, R1j, R2j, thetaj){

  Ni = nrow(Xi)

  Sigma1 <- create_sigma2(datlist = list(Zi, Zi), R1j)
  Sigma2 <- create_sigma2(datlist = list(Zi, Zi), R2j)

  beta1 <- c(beta1j[1, ], beta1j[2, ])
  beta2 <- c(beta2j[1, ], beta2j[2, ])

  Y <- c(Yi[, 1], Yi[, 2])

  Z <- direct.sum(Zi, Zi); X <- direct.sum(Xi, Xi)

  # draw bi over time (use posterior mean for now)
  b1ji = G1j %*% t(Z) %*% solve(Z %*% G1j %*% t(Z) + Sigma1) %*% (Y - X %*% beta1)
  b2ji = G2j %*% t(Z) %*% solve(Z %*% G2j %*% t(Z) + Sigma2) %*% (Y - X %*% beta2)

  mu1 = Xi %*% t(beta1j) + cbind(Zi[,1]*b1ji[1], Zi[,1]*b1ji[3]) + cbind(Zi[,2]*b1ji[2], Zi[,2]*b1ji[4])
  mu2 = Xi %*% t(beta2j) + cbind(Zi[,1]*b2ji[1], Zi[,1]*b2ji[3]) + cbind(Zi[,2]*b2ji[2], Zi[,2]*b2ji[4])

  lps = log(thetaj)

  lpdf1 = unlist(lapply(1:nrow(Xi), function(n) dmvnorm(Yi[n,], mean = mu1[n,], sigma = R1j, log = T) ))
  lpdf2 = unlist(lapply(1:nrow(Xi), function(n) dmvnorm(Yi[n,], mean = mu2[n,], sigma = R2j, log = T) ))

  logpyz = matrix(rep(lps, Ni), byrow=T, nrow=Ni) + cbind(cumsum(lpdf1), cumsum(lpdf2))

  logdenom = log(exp(logpyz[,1])+exp(logpyz[,2]))
  postz1 = exp(logpyz[,1] - logdenom)

  return(cbind(mu1, mu2, postz1))

}


calcpostzb0 = function(Xi, Zi, Yi, beta1j, beta2j, G1j, G2j, R1j, R2j, thetaj){

  Ni = nrow(Xi)


  Sigma1 <- create_sigma2_b0(datlist = list(Zi, Zi), R1j)
  Sigma2 <- create_sigma2_b0(datlist = list(Zi, Zi), R2j)

  beta1 <- c(beta1j[1, ], beta1j[2, ])
  beta2 <- c(beta2j[1, ], beta2j[2, ])

  Y <- c(Yi[, 1], Yi[, 2])

  Z <- direct.sum(Zi, Zi); X <- direct.sum(Xi, Xi)

  # draw bi over time (use posterior mean for now), BLUP
  b1ji = G1j %*% t(Z) %*% solve(Z %*% G1j %*% t(Z) + Sigma1) %*% (Y - X %*% beta1)
  b2ji = G2j %*% t(Z) %*% solve(Z %*% G2j %*% t(Z) + Sigma2) %*% (Y - X %*% beta2)

  mu1 = Xi %*% t(beta1j) + cbind(Zi[,1]*b1ji[1], Zi[,1]*b1ji[2]) #+ cbind(Zi[,2]*b1ji[2], Zi[,2]*b1ji[4])
  mu2 = Xi %*% t(beta2j) + cbind(Zi[,1]*b2ji[1], Zi[,1]*b2ji[2]) #+ cbind(Zi[,2]*b2ji[2], Zi[,2]*b2ji[4])

  lps = log(thetaj)

  lpdf1 = unlist(lapply(1:nrow(Xi), function(n) dmvnorm(Yi[n,], mean = mu1[n,], sigma = R1j, log = T) ))
  lpdf2 = unlist(lapply(1:nrow(Xi), function(n) dmvnorm(Yi[n,], mean = mu2[n,], sigma = R2j, log = T) ))

  logpyz = matrix(rep(lps, Ni), byrow=T, nrow=Ni) + cbind(cumsum(lpdf1), cumsum(lpdf2))

  logdenom = log(exp(logpyz[,1])+exp(logpyz[,2]))
  postz1 = exp(logpyz[,1] - logdenom)

  return(cbind(mu1, mu2, postz1))

}

####################################################################################
# shared X
calcpostzb0sX = function(t, Xi, xSi, Zi, Yi, beta1j, beta2j, betaSj, G1j, G2j, R1j, R2j, thetaj){

  beta1 <- c(beta1j[1, ], beta1j[2, ])
  beta2 <- c(beta2j[1, ], beta2j[2, ])
  beta_S <- c(betaSj[1, ], betaSj[2, ])

  # time-varying design matrices
  Zit <- matrix(Zi[1:t, ], nrow = t); Xit <- matrix(Xi[1:t, ], nrow = t)
  xSit <- matrix(xSi[1:t, ], nrow = t); Yit <- Yi[1:t, ]

  # var(Y)
  if(t == 1){
    Sigma1 = Sigma2 = R1j
  }else{
    Sigma1 = Sigma2 <- create_sigma2_b0(datlist = list(Zit, Zit), R1j)
  }

  # reformat design matrices and Y
  Z <- direct.sum(Zit[, 1], Zit[, 1]); X <- direct.sum(Xit, Xit); Xs <- direct.sum(xSit, xSit);
  Y <- rbind(Yit[, 1], Yit[, 2], use.names = F) %>% as.matrix

  # simulate b from b|Y for cluster 1
  b1ji = G1j %*% t(Z) %*% solve(Z %*% G1j %*% t(Z) + Sigma1) %*%
    (Y - (cbind(X, Xs) %*% c(beta1, beta_S)))

  varb1ji = G1j - G1j %*% t(Z) %*% solve(Z %*% G1j %*% t(Z) + Sigma1) %*% Z %*% G1j

  b1_sim <- rmvnorm(1, mean = b1ji, sigma = varb1ji)

  # simulate b from b|Y for cluster 2
  b2ji = G2j %*% t(Z) %*% solve(Z %*% G2j %*% t(Z) + Sigma2) %*%
    (Y - (cbind(X, Xs) %*% c(beta2, beta_S)))

  varb2ji = G2j - G2j %*% t(Z) %*% solve(Z %*% G2j %*% t(Z) + Sigma2) %*% Z %*% G2j

  b2_sim <- rmvnorm(1, mean = b2ji, sigma = varb2ji)

  # simulate random error
  eps1 = eps2 <- rmvnorm(1, mean = rep(0, t*2), sigma = Sigma1)

  Y1_sim = cbind(X, Xs) %*% c(beta1, beta_S) + Z %*% t(b1_sim) + t(eps1)
  Y2_sim = cbind(X, Xs) %*% c(beta2, beta_S) + Z %*% t(b2_sim) + t(eps2)

  c1fvc = Y1_sim[1:t]; c1dlco = Y1_sim[(t+1):(2*t)]
  c2fvc = Y2_sim[1:t]; c2dlco = Y2_sim[(t+1):(2*t)]

  return(cbind(c1fvc, c1dlco, c2fvc, c2dlco)[t,]) # vector of length 4

}

####################################################################################


eachperson = function(i, j){## for ith person jth draw in the one permuted large chain

  pd = which(refdat$pt_id == i)
  Yi = Yvec[pd,]
  Xi = xY[pd,]
  Zi = zY[pd,]

  beta1j = res$beta_1[j,,]
  beta2j = res$beta_2[j,,]
  G1j = res$G_1[j,,]
  G2j = res$G_2[j,,]
  R1j = res$R[j,,]
  R2j = R1j
  thetaj = res$theta[j,]

  mt = calcpostzb0bt(Xi, Zi, Yi, beta1j, beta2j, G1j, G2j, R1j, R2j, thetaj)
  return(mt)
}


eachperson_NP = function(i, l, j){## for ith person lth chain's jth draw

  pd = which(refdat$pt_id == i)
  Yi = Yvec[pd,]
  Xi = xY[pd,]
  if(ncol(zY)==1){
    Zi = matrix(zY[pd,], ncol = 1)
  } else {
    Zi = zY[pd,]
  }

  beta1j = matrix(res[j,l,grep("beta_1", dn$parameters)], nrow=2)
  beta2j = matrix(res[j,l,grep("beta_2", dn$parameters)], nrow=2)

  if(ncol(Zi) == 1){
    G1j = matrix(res[j,l,grep("G_1", dn$parameters)], ncol=2)
    G2j = matrix(res[j,l,grep("G_2", dn$parameters)], ncol=2)
  } else {
    G1j = matrix(res[j,l,grep("G_1", dn$parameters)], ncol=4)
    G2j = matrix(res[j,l,grep("G_2", dn$parameters)], ncol=4)
  }

  R1j = matrix(res[j,l,grep("R", dn$parameters)], ncol=2)
  R2j = R1j
  thetaj = res[j,l,grep("theta", dn$parameters)]

  if(ncol(Zi) == 1){
    mt = calcpostzb0(Xi, Zi, Yi, beta1j, beta2j, G1j, G2j, R1j, R2j, thetaj)
  } else {
    mt = calcpostzb0bt(Xi, Zi, Yi, beta1j, beta2j, G1j, G2j, R1j, R2j, thetaj)
  }

  return(mt)
}

lfun = function(x){
  s = sort(x)
  n = length(x)
  return(s[ceiling(n*0.025)])
}
ufun = function(x){
  s = sort(x)
  n = length(x)
  return(s[ceiling(n*0.975)])
}


####################################################################################
eachperson_NP_xS = function(i, l, j){## for ith person lth chain's jth draw
  # shared X, random intercept only

  pd = which(refdat$pt_id == i)
  Yi = Yvec[pd, ]
  Xi = xY[pd, ]
  xSi = xS[pd, ]
  Zi = zY[pd,]

  beta1j = matrix(res[j,l,grep("beta_1", dn$parameters)], nrow=2)
  beta2j = matrix(res[j,l,grep("beta_2", dn$parameters)], nrow=2)

  betaSj = matrix(res[j,l,grep("betaS", dn$parameters)], nrow=2)

  G1j = matrix(res[j,l,grep("G_1", dn$parameters)], ncol=2)
  G2j = matrix(res[j,l,grep("G_2", dn$parameters)], ncol=2)

  R1j = matrix(res[j,l,grep("R", dn$parameters)], ncol=2)
  R2j = R1j

  mt = lapply(1:nrow(Yi), function(t) calcpostzb0sX(t, Xi, xSi, Zi, Yi, beta1j, beta2j, betaSj, G1j, G2j, R1j, R2j))
  mt = t(simplify2array(mt)) # Ti x 4

  return(mt)
}


####################################################################################



drawepsilon_NP = function(xY, zY, Yvec, l, j){## for ith person lth chain's jth draw


  beta1j = matrix(res[j,l,grep("beta_1", dn$parameters)], nrow=2)
  beta2j = matrix(res[j,l,grep("beta_2", dn$parameters)], nrow=2)
  b1j = matrix(res[j,l,grep("b_1", dn$parameters)], ncol=4)
  b2j = matrix(res[j,l,grep("b_2", dn$parameters)], ncol=4)
  b1j = b1j[refdat$ID,]
  b2j = b2j[refdat$ID,]
  mu1 = xY %*% t(beta1j) + cbind(zY[,1]*b1j[,1], zY[,1]*b1j[,3]) + cbind(zY[,2]*b1j[,2], zY[,2]*b1j[,4])
  mu2 = xY %*% t(beta2j) + cbind(zY[,1]*b2j[,1], zY[,1]*b2j[,3]) + cbind(zY[,2]*b2j[,2], zY[,2]*b2j[,4])

  epsilon1 = Yvec - mu1
  epsilon2 = Yvec - mu2
  return(list(epsilon1, epsilon2))
}

### 2023 Jan 9th: real time update with functional cluster membership
# Jan 10th edits: b1i0 vector into columns

library(mvtnorm) # dmvnorm, rmvnorm
library(boot) # inv.logit
gij_fun = function(Yalli, xYi, xSi, beta1j, beta2j, betaSj, G1j, G2j, Rj, Zi, alphaj, b1, b2){
  # in the arguments, i index person, j index posterior draw

  # number of visits
  ni = nrow(Yalli)

  # cluster 1
  mu1 = xYi %*% beta1j + xSi %*% betaSj + matrix(rep(b1,ni), byrow=T, nrow=ni) # matrix of ni x 2, columns correspond to FVC and DLCO
  pY1_individual = unlist(lapply(1:ni, function(t) dmvnorm(x = Yalli[t,], mean = mu1[t,], sigma = Rj))) # vector of length ni
  #pY1 = cumprod(pY1_individual)

  # cluster 2
  mu2 = xYi %*% beta2j + xSi %*% betaSj + matrix(rep(b2,ni), byrow=T, nrow=ni)
  pY2_individual = unlist(lapply(1:ni, function(t) dmvnorm(x = Yalli[t,], mean = mu2[t,], sigma = Rj))) # vector of length ni
  #pY2 = cumprod(pY2_individual)

  # P(cluster 1)
  p1 = inv.logit(as.matrix(Zi) %*% matrix(alphaj, ncol = 1)) # vector of length ni

  # git jth posterior
  #gij = pY1 * p1 / ( pY1 * p1 +  pY2 * (1 - p1)) # vector of length ni
  gij = 1 / ( 1 +  cumprod(pY2_individual/pY1_individual) * (1/p1 - 1))
  return(gij)
}

# calculate likelihood for cluster 1 & 2
L_mvn <- function(t, Yalli, xYi, xSi, beta1j, beta2j, betaSj, G1j, G2j, Rj){

  In2 <- diag(2)

  zmat <- In2 %x% matrix(rep(1, t), byrow = T, nrow = t) # matrix of t x 2, columns correspond to FVC and DLCO
  yvar_c1 <- zmat %*% G1j %*% t(zmat) + Rj %x% diag(t)
  yvar_c2 <- zmat %*% G2j %*% t(zmat) + Rj %x% diag(t)

  YALLi <- rbind(Yalli[1:t, 1], Yalli[1:t, 2], use.names = F) %>% as.matrix

  # cluster 1
  mu1 = xYi %*% beta1j + xSi %*% betaSj
  Mu1 <- matrix(c(mu1[1:t, 1], mu1[1:t, 2]), ncol = 1)

  # cluster 2
  mu2 = xYi %*% beta2j + xSi %*% betaSj
  Mu2 <- matrix(c(mu2[1:t, 1], mu2[1:t, 2]), ncol = 1)

  c1dmvn <- dmvnorm(x = c(YALLi), mean = c(Mu1), sigma = yvar_c1)
  c2dmvn <- dmvnorm(x = c(YALLi), mean = c(Mu2), sigma = yvar_c2)

  return(list(c1dmvn, c2dmvn))

}

# gij calculated using multivariate normal likelihood
gij_mvn = function(Yalli, xYi, xSi, beta1j, beta2j, betaSj, G1j, G2j, Rj, Zi, alphaj){
  # in the arguments, i index person, j index posterior draw

  # number of visits
  ni = nrow(Yalli)

  # individual's likelihood
  pY_individual = t(matrix(unlist(lapply(1:ni, function(t) L_mvn(t, Yalli, xYi, xSi, beta1j, beta2j, betaSj, G1j, G2j, Rj))), nrow = 2))

  # log likelihood
  L_pY1 <- log(pY_individual[, 1]); L_pY2 <- log(pY_individual[, 2])

  # P(cluster 1)
  p1 = inv.logit(as.matrix(Zi) %*% matrix(alphaj, ncol = 1)) # vector of length ni

  # git jth posterior
  #gij = pY1 * p1 / ( pY1 * p1 +  pY2 * (1 - p1)) # vector of length ni
  gij = 1 / ( 1 +  exp(L_pY2 - L_pY1) * (1/p1 - 1))

  return(gij)

}


if(0){
  i = 2
  id = idvec[i]

  pidind <- which(refdat$pt_id == 20)
  Yalli = refdat[pidind, list(fvc,dlco)]
  xYi <- xY[pidind, ]
  xSi <- xS[pidind, ]
  Zi <- Z[pidind, ]


  p1 <- inv.logit(as.matrix(Zi) %*% matrix(alphaj, ncol = 1))

  #tmp = cbind(Yalli,p1,mu1,mu2)
  #colnames(tmp) = c("FVC","DLCO","theta","muFVC1","muDLCO1","muFVC2","muDLCO2")

  beta1j <- t(matrix(apply(beta1, 2, mean), nrow = 2))
  beta2j <- t(matrix(apply(beta2, 2, mean), nrow = 2))
  betaSj <- t(matrix(apply(betaS, 2, mean), nrow = 2))

  G1j <- matrix(apply(G1, 2, mean), nrow = 2)
  G2j <- matrix(apply(G2, 2, mean), nrow = 2)
  Rj <- matrix(apply(R, 2, mean), nrow = 2)
  alphaj <- apply(alpha, 2, mean)
  gij_mvn(Yalli, xYi, xSi, beta1j, beta2j, betaSj, G1j, G2j, Rj, Zi, alphaj)

}

# 2023/06/24

####################################################################################
# shared X sampling b|Y(1: t-1) for Yt
calcpostzb0sX_pred = function(t, Xi, xSi, Zi, Yi, beta1j, beta2j, betaSj, G1j, G2j, R1j, R2j, thetaj){
  
  beta1 <- c(beta1j[1, ], beta1j[2, ])
  beta2 <- c(beta2j[1, ], beta2j[2, ])
  beta_S <- c(betaSj[1, ], betaSj[2, ])
  
  # time-varying design matrices
  Zit <- matrix(Zi[1:t, ], nrow = t); Xit <- matrix(Xi[1:t, ], nrow = t)
  xSit <- matrix(xSi[1:t, ], nrow = t); Yit <- Yi[1:t, ]
  
  # var(Y)
  if(t == 1){
    Sigma1 = Sigma2 = R1j
  }else{
    Sigma1 = Sigma2 <- create_sigma2_b0(datlist = list(Zit, Zit), R1j)
  }
  
  # reformat design matrices and Y
  Z <- direct.sum(Zit[, 1], Zit[, 1]); X <- direct.sum(Xit, Xit); Xs <- direct.sum(xSit, xSit);
  Y <- rbind(Yit[, 1], Yit[, 2], use.names = F) %>% as.matrix
  
  if(t == 1){
    b1_sim <- rmvnorm(1, mean = rep(0,nrow(G1j)), sigma = G1j)
    b2_sim <- rmvnorm(1, mean = rep(0,nrow(G2j)), sigma = G2j)
  } else {
    Zit_prev <- matrix(Zi[1:(t-1), ], nrow = t-1); Xit_prev <- matrix(Xi[1:(t-1), ], nrow = t-1)
    xSit_prev <- matrix(xSi[1:(t-1), ], nrow = t-1); Yit_prev <- Yi[1:(t-1), ]
    Sigma1_prev = Sigma2_prev <- create_sigma2_b0(datlist = list(Zit_prev, Zit_prev), R1j)
    Z_prev <- direct.sum(Zit_prev[, 1], Zit_prev[, 1]); X_prev <- direct.sum(Xit_prev, Xit_prev); Xs_prev <- direct.sum(xSit_prev, xSit_prev);
    Y_prev <- rbind(Yit_prev[, 1], Yit_prev[, 2], use.names = F) %>% as.matrix
    
    
    # simulate b from b|Y for cluster 1
    b1ji = G1j %*% t(Z_prev) %*% solve(Z_prev %*% G1j %*% t(Z_prev) + Sigma1_prev) %*%
      (Y_prev - (cbind(X_prev, Xs_prev) %*% c(beta1, beta_S)))
    
    varb1ji = G1j - G1j %*% t(Z_prev) %*% solve(Z_prev %*% G1j %*% t(Z_prev) + Sigma1_prev) %*% Z_prev %*% G1j
    
    b1_sim <- rmvnorm(1, mean = b1ji, sigma = varb1ji)
    
    # simulate b from b|Y for cluster 2
    b2ji = G2j %*% t(Z_prev) %*% solve(Z_prev %*% G2j %*% t(Z_prev) + Sigma2_prev) %*%
      (Y_prev - (cbind(X_prev, Xs_prev) %*% c(beta2, beta_S)))
    
    varb2ji = G2j - G2j %*% t(Z_prev) %*% solve(Z_prev %*% G2j %*% t(Z_prev) + Sigma2_prev) %*% Z_prev %*% G2j
    
    b2_sim <- rmvnorm(1, mean = b2ji, sigma = varb2ji)
  }
  
  # simulate random error
  eps1 = eps2 <- rmvnorm(1, mean = rep(0, t*2), sigma = Sigma1)
  
  Y1_sim = cbind(X, Xs) %*% c(beta1, beta_S) + Z %*% t(b1_sim) + t(eps1)
  Y2_sim = cbind(X, Xs) %*% c(beta2, beta_S) + Z %*% t(b2_sim) + t(eps2)
  
  c1fvc = Y1_sim[1:t]; c1dlco = Y1_sim[(t+1):(2*t)]
  c2fvc = Y2_sim[1:t]; c2dlco = Y2_sim[(t+1):(2*t)]
  
  return(cbind(c1fvc, c1dlco, c2fvc, c2dlco)[t,]) # vector of length 4
  
}

####################################################################################
# sampling b|Y(1: t-1) for Yt
eachperson_NP_xS_pred = function(i, l, j){## for ith person lth chain's jth draw
  # shared X, random intercept only
  
  pd = which(refdat$pt_id == i)
  Yi = Yvec[pd, ]
  Xi = xY[pd, ]
  xSi = xS[pd, ]
  Zi = zY[pd,]
  
  beta1j = matrix(res[j,l,grep("beta_1", dn$parameters)], nrow=2)
  beta2j = matrix(res[j,l,grep("beta_2", dn$parameters)], nrow=2)
  
  betaSj = matrix(res[j,l,grep("betaS", dn$parameters)], nrow=2)
  
  G1j = matrix(res[j,l,grep("G_1", dn$parameters)], ncol=2)
  G2j = matrix(res[j,l,grep("G_2", dn$parameters)], ncol=2)
  
  R1j = matrix(res[j,l,grep("R", dn$parameters)], ncol=2)
  R2j = R1j
  
  mt = lapply(1:nrow(Yi), function(t) calcpostzb0sX_pred(t, Xi, xSi, Zi, Yi, beta1j, beta2j, betaSj, G1j, G2j, R1j, R2j))
  mt = t(simplify2array(mt)) # Ti x 4
  
  return(mt)
}




















