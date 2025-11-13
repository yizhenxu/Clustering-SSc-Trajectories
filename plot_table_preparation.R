if(0){
  load("S:/TrajectoryAnalysis/STAN/res/mod3res_shrR_nR_b0_221115.RData")
  load("S:/TrajectoryAnalysis/STAN/res/mod2res_shrR_nR_b0_221115.RData")
  load("S:/TrajectoryAnalysis/STAN/res/mod2res_shrR_nR_b0_221024.RData")
  load("S:/TrajectoryAnalysis/STAN/res/mod1res_shrR_nR_b0_221115.RData")
  load("S:/TrajectoryAnalysis/STAN/res/mod2res_shrR_nR_b0_221117_rnapol.RData")
  load("S:/TrajectoryAnalysis/STAN/res/mod2res_shrR_nR_b0_221117_scl70.RData")
  load("S:/TrajectoryAnalysis/STAN/res/mod2res_shrR_nR_b0_221127_all.RData")
  load("S:/TrajectoryAnalysis/STAN/res/mod2res_sharedBeta_221129.RData")
  load("S:/TrajectoryAnalysis/STAN/res/mod2res_sharedBeta_221213_scl70.RData")
  load("S:/TrajectoryAnalysis/STAN/res/mod2res_sharedBeta_230813_scl70.RData")
}

load("S:/TrajectoryAnalysis/STAN/res/mod2res_sharedBeta_230813_scl70_functionalTheta.RData")

dn = dimnames(res)

lp = res[, , which(dn$parameters == "lp__")]
odch = order(apply(lp, 2, max), decreasing = T)

chnm = odch[1]

betaS = res[,chnm, grep("betaS",dn$parameters )]; dimnames(betaS)$parameters

beta1 = res[,chnm, grep("beta_1",dn$parameters )]; dimnames(beta1)$parameters

beta2 = res[,chnm, grep("beta_2",dn$parameters )]; dimnames(beta2)$parameters

#beta3 = res[,chnm, grep("beta_3",dn$parameters )]; dimnames(beta3)$parameters

G1 = res[,chnm,grep("G_1", dn$parameters)]; dimnames(G1)$parameters

G2 = res[,chnm,grep("G_2", dn$parameters)]; dimnames(G2)$parameters

#G3 = res[,chnm,grep("G_3", dn$parameters)]; dimnames(G3)$parameters

R = res[,chnm,grep("R", dn$parameters)]; dimnames(R)$parameters

b1 = res[,chnm, grep("b_1",dn$parameters )]; dimnames(b1)$parameters

b2 = res[,chnm, grep("b_2",dn$parameters )]; dimnames(b2)$parameters

#b3 = res[,chnm, grep("b_3",dn$parameters )]; dimnames(b3)$parameters

alpha = res[,chnm, grep("alpha",dn$parameters )];

summ_report = function(x){
  m = mean(x)
  l = length(x)
  ci = sort(x)[round(l*c(0.025,0.975))]
  r = c(m, ci)
  r = round(r,2);
  return(paste0(r[1], " (", r[2], ", ",r[3],")" ))

}

summ_original = function(x){ #
  m = mean(x)
  l = length(x)
  ci = sort(x)[round(l*c(0.025,0.975))]
  return(c(m, ci))
}

summ_standardized = function(x){ # standardized
  m = mean(x)
  s = sd(x)
  l = length(x)
  ci = sort(x)[round(l*c(0.025,0.975))]
  return(c(m/s, ci/s))
}

summ_byYsd = function(x,s){ # standardized
  m = mean(x)
  #s = sd(x)
  l = length(x)
  ci = sort(x)[round(l*c(0.025,0.975))]
  return(c(m/s, ci/s))
}

# summaries of betas
# qqnorm for b1 and b2
# summaries of Gs and the shared R by cluster

# summaries of X by cluster membership
# population trajectories by covariate subgroup and cluster
# cluster probabilities using all information
# spaghetties by cluster

### check psrf
if(0){
  #install.packages("stableGR")
  library(stableGR)
 
  betalist = cbind(beta1, beta2, betaS)
  Glist = cbind(G1[,c(1,2,4)], G2[,c(1,2,4)])
  
  stable.GR(list(betalist))
  stable.GR(list(Glist))
  stable.GR(list(R[,c(1,2,4)]))
  stable.GR(list(alpha))
}
