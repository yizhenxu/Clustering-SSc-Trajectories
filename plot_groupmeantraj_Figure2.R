#### plot mean trajectory by group
# adapt from plot_meantraj_by_covariate_sharedBeta.R

setwd("S://TrajectoryAnalysis//STAN")

source("./BMC submission 2024 March/plot_table_preparation.R") # related to task 2

# load functions: git_fun
source("Functions.R")


nit = 5000
nwup = 2000
seednum = 1
nc = 10 # number of chains

min_Ni = 2 # at least min_Ni obs for each person
subsetN = F # take a subset of patients


library(plyr)
library(dplyr)
library(codetools)
library(rstan)
library(splines)
library(data.table)

source("dataprocessing3_0-10_scl70.R")

meanage <- mean(refdat[!duplicated(refdat$pt_id), ]$ageonset, na.rm = T)

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


# Z alpha - defining matrix for function theta
#Zmat <- refdat %>% filter(!duplicated(pt_id)) %>% select(all_of(factorname))
Zmat <- refdat %>% select(all_of(factorname))

Z = cbind(1, Zmat)
colnames(Z)[1] <- "Intercept"

# second column needed to create sigma matrix
zY = cbind(1, refdat$YearsSinceOnset)

### up to this point -- preparation

# beta

# plot cluster specific trajectories

# plot cluster specific trajectories
CIbyrow <- function(data){apply(data, 1, quantile, c(0.025, 0.975))}

k1ind <- grep("\\[1", colnames(beta1))

beta1.pm <- apply(beta1, 2, mean)
beta2.pm <- apply(beta2, 2, mean)

beta1.pm.fvc <- beta1.pm[k1ind]; beta1.pm.dlco <- beta1.pm[-k1ind]
beta2.pm.fvc <- beta2.pm[k1ind]; beta2.pm.dlco <- beta2.pm[-k1ind]

betaS.pm <- apply(betaS, 2, mean)

k1indS <- grep("\\[1", names(betaS.pm))
betaS.pm.fvc <- betaS.pm[k1indS]; betaS.pm.dlco <- betaS.pm[-k1indS]

# design matrix

#nsres <- ns(refdat$YearsSinceOnset, df = 3, Boundary.knots = c(0, 20))
#knots <- attr(nsres, "knots")

yearvec <- seq(0, 10, by = 0.1)
tmat <- model.matrix(~ ns(yearvec, knots =5, Boundary.knots = c(0, 10)))


# function to plot mean trajectories

factorname <- c("male", "AArace","diffuse", "late_ageonset")

plot_groupmeantraj <- function(xvec){

  #ptitle <- paste0(paste(factorname[-3][as.logical(xvec[-3])], collapse = ", "),
  #                 paste(", onset age", round(xvec[3])))
  #xvec[3] <- xvec[3] - meanage #mean_onsetage
  ptitle <- paste0(factorname, xvec , collapse = ", ")

  xmat <- as.data.table(matrix(xvec, nrow = 1))
  colnames(xmat) <- factorname

  xs = c()
  xsname = c()
  for(i in 1:ncol(xmat)){
    for(j in 1:ncol(tmat)){
      xsname = c(xsname, paste0(factorname[i],"_", (j-1)))
      xs = cbind(xs, as.numeric(xmat[, get(factorname[i])])*tmat[,j])
    }
  }

  colnames(xs) = xsname

  
  # cluster specific X
  xy = tmat
  colnames(xy) <- c("Intercept", paste0("YearsSinceOnset",1:2))
  
  # plot trajectories
  xn <- length(yearvec)

  
  k1ind <- grep("\\[1", colnames(beta1))
  k1indS <- grep("\\[1", colnames(betaS))
  fvcmat1 <- xy %*% t(beta1[,k1ind]) + xs %*% t(betaS[,k1indS])
  fvcmat2 <- xy %*% t(beta2[,k1ind]) + xs %*% t(betaS[,k1indS])
  dlcomat1 <- xy %*% t(beta1[,-k1ind]) + xs %*% t(betaS[,-k1indS])
  dlcomat2 <- xy %*% t(beta2[,-k1ind]) + xs %*% t(betaS[,-k1indS])
  
  # plot cluster specific trajectories
  CIs = cbind(CIbyrow(fvcmat1), CIbyrow(dlcomat1), CIbyrow(fvcmat2), CIbyrow(dlcomat2))
  
  pdat <- data.frame(group = c(rep("Fast progressor", 2*xn), rep("Stable", 2*xn)),
                     Measure = rep(c(rep("pFVC", xn), rep("pDLCO", xn)), 2),
                     year = rep(yearvec, 4),
                     xbeta = c(xy %*% beta1.pm.fvc + xs %*% betaS.pm.fvc, xy %*% beta1.pm.dlco + xs %*% betaS.pm.dlco,
                               xy %*% beta2.pm.fvc + xs %*% betaS.pm.fvc, xy %*% beta2.pm.dlco + xs %*% betaS.pm.dlco),
                     lower = CIs[1,],
                     upper = CIs[2,])

  library(ggplot2)
  ggplot(pdat, aes(x = year, group = Measure, fill = Measure)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    geom_line(aes(y = xbeta, color = Measure)) + #theme(text = element_text(size = 20))+
    facet_wrap(~group) + labs(x = "Years since onset", y = "Estimated trajectory")#, title = ptitle


}


factorname
xvec <- c(0,0,0,0)
#xvec <- c(0, 0, 24, 1, 0, 1, 0)

pdf("./BMC submission 2024 March/Figure 2.pdf", width = 7, height = 3.5)
plot_groupmeantraj(xvec)
dev.off()


