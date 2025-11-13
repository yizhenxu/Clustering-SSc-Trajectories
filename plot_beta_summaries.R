### Nightingale of fix effect betas by cluster and outcome

################################################################################
### load estimated parameters
setwd("S:/TrajectoryAnalysis/STAN")

source("./BMC submission 2024 March/plot_table_preparation.R")

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
source("S:\\TrajectoryAnalysis\\STAN\\geom_effect.R")


### Nightingale plot
setDT(plotd)
p1 <- ggplot(data = plotd[Y=="FVC" ,], aes(x = beta, y = reorder(Xname,Xorder))) +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.5, alpha=0.6)+
  geom_effect(
    ggplot2::aes(
      xmin = bmin,
      xmax = bmax,
      colour =as.factor( Cluster )
    ),
    position = ggstance::position_dodgev(height = 0.5)
  )+
  labs(title="FVC")+theme(legend.position="bottom")

p2 <- ggplot(data = plotd[Y=="DLCO" ,], aes(x = beta, y = reorder(Xname,Xorder))) +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.5, alpha=0.6)+
  geom_effect(
    ggplot2::aes(
      xmin = bmin,
      xmax = bmax,
      colour =as.factor( Cluster )
    ),
    position = ggstance::position_dodgev(height = 0.5)
  )+
  labs(title="DLCO")+theme(legend.position="bottom")

png("./Results221004/plot_betas_by_Cluster_2_230813.png", width = 700, height = 500)
plot(arrangeGrob(p1 ,
                 p2 ,
                 nrow=1))
dev.off()


