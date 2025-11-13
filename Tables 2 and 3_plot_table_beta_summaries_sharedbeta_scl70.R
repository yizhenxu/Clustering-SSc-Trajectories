### Nightingale of fix effect betas by cluster and outcome

################################################################################
### load estimated parameters
setwd("S://TrajectoryAnalysis//STAN")

# notes for Yizhen: script with new function gij_mvn

# load results: beta1, beta2, betaS, G1, G2, R
source("./BMC submission 2024 March/plot_table_preparation.R") # related to task 2

# load functions: git_fun
source("Functions.R")

#setwd("./workspace/Storage/yxu143/persistent/Scleroderma clustering")
#load("mod2res_sharedBeta_221213_scl70.RData")
#load("S:/TrajectoryAnalysis/STAN/res/mod2res_sharedBeta_221213_scl70_functionalTheta.RData")

#load("S:/TrajectoryAnalysis/STAN/res/mod2res_sharedBeta_230813_scl70_functionalTheta.RData")

################################################################################
### load design matrix xY

library(plyr)
library(dplyr)
library(codetools)
library(rstan)
library(splines)
library(data.table)

nit = 5000
nwup = 2000
seednum = 1
nc = 10 # number of chains

min_Ni = 2 # at least min_Ni obs for each person
subsetN = F # take a subset of patients

# if subsetN == T, specify the seed and number of ppl to subset
#subset_seed = 123; Nsub = 500

# data processing script
source("dataprocessing3_0-10_scl70.R")


tmp = refdat[, sum(is.na(diffuse) + is.na(scl70)), by = pt_id]
length(unique(tmp$pt_id[tmp$V1 == 0])) # 284

modIDset = tmp$pt_id[tmp$V1 == 0]
refdat = refdat[pt_id %in% modIDset,] # 1046 ppl

# define lsub for stan
ltab <- refdat %>% group_by(pt_id) %>% summarize(lsub = length(pt_id))
lsub <- ltab$lsub

pidtab <- table(refdat$pt_id)
nid <- length(pidtab)

refdat$ID <- rep(1:nid, pidtab)

# define
Yvec <- refdat[, c("fvc", "dlco")]

Tmat <- model.matrix(~ ns(YearsSinceOnset, knots = 5, Boundary.knots = c(0, 10)), refdat)

factorname <- c("male", "AArace","diffuse", "late_ageonset")
#c("male", "AArace", "ageonset", "diffuse", "aca", "rnapol", "scl70")
Xmat <- refdat %>% select(all_of(factorname))

# shared X
xSname = c(); xS = c()
for(i in 1:ncol(Xmat)){
  for(j in 1:ncol(Tmat)){
    xSname = c(xSname, paste0(factorname[i],"_", (j-1)))
    xS = cbind(xS, as.numeric(Xmat[,get(factorname[i])])*Tmat[,j])
  }
}
colnames(xS) = xSname

# cluster specific X
xY = Tmat
colnames(xY) <- c("Intercept", paste0("YearsSinceOnset",1:2))

################################################################################
################################################################################
### create table for plot for 2 clusters

beta1m = apply(beta1, 2, summ_original)
beta2m = apply(beta2, 2, summ_original)

yn = c("FVC", "DLCO")
L = 2
plotd = matrix(NA, ncol = 7, nrow = ncol(xY)*length(yn)*L )
colnames(plotd) = c("Xname", "Xorder", "Y", "beta", "bmin","bmax","Cluster")
plotd = as.data.frame(plotd)
plotd[,4:6] = rbind(t(beta1m), t(beta2m))

betaelem_X = matrix(rep(colnames(xY),2) , byrow=T, nrow=2)
betaelem_Xorder =  matrix(rep(ncol(xY):1,2) , byrow=T, nrow=2)
betaelem_Y = rbind(rep(yn[1], ncol(xY)), rep(yn[2], ncol(xY)))

plotd$Xname = rep(c(betaelem_X), L)
plotd$Xorder =  rep(c(betaelem_Xorder), L)
plotd$Y = rep(c(betaelem_Y), L)

plotd$Cluster = rep(1:2, each = ncol(xY)*length(yn) )

library(ggplot2)
library(magrittr)
library(gridExtra)
library(grid)
#install.packages("ggstance")
library(ggstance)
#source("S:\\TrajectoryAnalysis\\STAN\\geom_effect.R")


source("./geom_effect.R")

plotd$Xname[plotd$Xname == "YearsSinceOnset1"] = "Spline 1" #"Years Since Onset 1"
plotd$Xname[plotd$Xname == "YearsSinceOnset2"] = "Spline 2" #"Years Since Onset 2"
plotd$Y[plotd$Y == "FVC"] = "pFVC" #"Years Since Onset 1"
plotd$Y[plotd$Y == "DLCO"] = "pDLCO" #"Years Since Onset 2"

plotd$Cluster[plotd$Cluster == "1"] = "Fast progressor"
plotd$Cluster[plotd$Cluster == "2"] = "Stable"
plotd$Cluster = as.factor( plotd$Cluster )

### Nightingale plot
setDT(plotd)
library(ggstance)

p_beta <- ggplot(data = plotd, aes(x = beta, y = reorder(Xname, Xorder))) +
  geom_vline(xintercept = 0, linetype = "solid", size = 0.5, alpha=0.6)+
  geom_effect(
    ggplot2::aes(
      xmin = bmin,
      xmax = bmax,
      colour = Cluster),
    position = ggstance::position_dodgev(height = 0.5))+
  labs(title="", x = "", y = "") +
  theme(legend.position="bottom", legend.title=element_blank()) +
  facet_wrap(~Y)

pdf("./BMC submission 2024 March/plot_betas_by_Cluster_2.pdf", width = 7, height = 3)
p_beta
dev.off()

################################################################################
### create table for plot for shared parameters

beta1m = apply(betaS, 2, summ_original)

yn = c("FVC", "DLCO")
L=1
plotd = matrix(NA, ncol = 7, nrow = ncol(xS)*length(yn) )
colnames(plotd) = c("Xname", "Xorder", "Y", "beta", "bmin","bmax","Cluster")
plotd = as.data.frame(plotd)
plotd[,4:6] = t(beta1m)

cn0 = c("Male", "Black", "Diffuse", "Late Onset")
cn = c()
for(i in cn0){
  cn = c(cn, paste0(i, " ", c("", " x Spline 1", " x Spline 2")))
}

betaelem_X = matrix(rep(cn,2) , byrow=T, nrow=2)#matrix(rep(colnames(xS),2) , byrow=T, nrow=2)
betaelem_Xorder =  matrix(rep(ncol(xS):1,2) , byrow=T, nrow=2)
betaelem_Y = rbind(rep(yn[1], ncol(xS)), rep(yn[2], ncol(xY)))

plotd$Xname = c(betaelem_X)
plotd$Xorder =  c(betaelem_Xorder)
plotd$Y = c(betaelem_Y)

library(ggplot2)
library(magrittr)
library(gridExtra)
library(grid)
install.packages("ggstance")
library(ggstance)
source("S:\\TrajectoryAnalysis\\STAN\\geom_effect.R")

source("./workspace/Storage/yxu143/persistent/Scleroderma/geom_effect.R")

plotd$Y[plotd$Y == "FVC"] = "pFVC" #"Years Since Onset 1"
plotd$Y[plotd$Y == "DLCO"] = "pDLCO" #"Years Since Onset 2"

### Nightingale plot
setDT(plotd)
p_beta_c <- ggplot(data = plotd, aes(x = beta, y = reorder(Xname,Xorder))) +
  geom_vline(xintercept = 0, linetype = "solid", size = 0.5, alpha=0.6)+
  geom_effect(
    ggplot2::aes(
      xmin = bmin,
      xmax = bmax
    ),
    position = ggstance::position_dodgev(height = 0.5)
  )+
  labs(title="", x = "", y = "")+theme(legend.position="bottom") +
  facet_wrap(~Y)


pdf("./BMC submission 2024 March/plot_betas_by_Cluster_2_shared.pdf", width = 7, height = 4)
p_beta_c
dev.off()


##########################################################
##########################################################
##########################################################
##########################################################
### CREATE TABLES

beta1m = apply(beta1, 2, summ_original)
beta2m = apply(beta2, 2, summ_original)

tab = cbind( c(c("1", "", ""), c("2", "", "")),
             rep(c("Intercept", "Years Since Onset 1","Years Since Onset 2"), 2),
             round(matrix(c(beta1m, beta2m), ncol = 6, byrow = T), 2) )
colnames(tab) = c("Cluster", "Variable", paste0("FVC_", c("Mean", "L", "U")), paste0("DLCO_", c("Mean", "L", "U")))
noquote(tab)


beta1m = apply(betaS, 2, summ_original)

cn0 = c("Male", "African American", "Diffuse", "Late Onset")
cn = c()
for(i in cn0){
  cn = c(cn, paste0(i, " ", 0:2))
}

tab =  rbind(tab, cbind(c("Shared", rep("", length(cn)-1) ), cn, matrix( round(beta1m,2), byrow = T, ncol = 6) ) )
noquote(tab)


save(tab, file = "table_beta_alpha.RData")


summ_mat = function(G1, ndim, varname, frontword){
  summ = apply(G1, 2, summ_original)
  r = apply(summ, 2, function(x) {x = round(x,2); return(paste0(x[1], " (", x[2], ", ", x[3],")" ))})
  m = matrix(r, ncol = ndim)
  colnames(m) = rownames(m) = varname
  print(noquote(frontword))
  return(noquote(m))
}


G1m = summ_mat(G1, 2, varname=c("FVC b0", "DLCO b0"), frontword = "Cluster 1 G:")
G2m = summ_mat(G2, 2, varname=c("FVC b0", "DLCO b0"), frontword = "Cluster 2 G:")
Rm = summ_mat(R, 2, varname=c("FVC", "DLCO"), frontword = "Shared R:")



alpham = apply(alpha, 2, summ_original)
tab_alpha = t(round(alpham, 2))
rownames(tab_alpha) = c("Intercept", "Male", "African American", "Diffuse", "Late Onset")

tab_alphaGR = matrix("", ncol = 3, nrow = 14)
tab_alphaGR[1:3,3] = G1m[c(1,2,4)]
tab_alphaGR[4:6,3] = G2m[c(1,2,4)]
tab_alphaGR[7:9,3] = Rm[c(1,2,4)]
tab_alphaGR[10:14,3] = paste0(tab_alpha[,1], ' (', tab_alpha[,2], ", ", tab_alpha[,3], ")")
tab_alphaGR[,2] = c(rep(c("G11", "G12", "G22"),2),c("R11", "R12", "R22"),  c("Intercept", "Male", "African American", "Diffuse", "Late Onset"))
tab_alphaGR[,1] = c("Fast progressor", rep("", 2),
                    "Stable", rep("", 2),
                    "Shared R", rep("", 2),
                    "Latent Coefficients", rep("", 4))
colnames(tab_alphaGR) = c("Cluster", "Variable", "Mean (95% CI)")
noquote(tab_alphaGR)

save(tab_alphaGR, file = "table_alphaGR.RData")

