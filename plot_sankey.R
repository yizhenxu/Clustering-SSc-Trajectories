
setwd("S://TrajectoryAnalysis//STAN")

# notes for Yizhen: script with new function gij_mvn

# load results: beta1, beta2, betaS, G1, G2, R
source("plot_table_preparation.R") # related to task 2

# load functions: git_fun
source("Functions.R")

### Prep

library(dplyr)
library(codetools)
library(rstan)
library(splines)
library(data.table)
library(mvtnorm)

nit = 5000
nwup = 2000
seednum = 1
nc = 10 # number of chains

min_Ni = 2 # at least min_Ni obs for each person
subsetN = F # take a subset of patients

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


# Z alpha - defining matrix for function theta
#Zmat <- refdat %>% filter(!duplicated(pt_id)) %>% select(all_of(factorname))
Zmat <- refdat %>% select(all_of(factorname))

Z = cbind(1, Zmat)
colnames(Z)[1] <- "Intercept"

# second column needed to create sigma matrix
zY = cbind(1, refdat$YearsSinceOnset)

### up to this point -- preparation

###########################################################################
###########################################################################
# calculating theta = P(cluster 1)

Ztmp1 <- refdat %>% filter(!duplicated(pt_id)) %>% select(all_of(factorname))
Ztmp2 = cbind(1, Ztmp1)

# matrix of N x npost  (289 x 3000)
library(boot) # inv.logit
alltheta = inv.logit(as.matrix(Ztmp2) %*% t(alpha))
#hist(alltheta)

###########################################################################
idvec = unique(refdat$pt_id)
#load("git_ALL_mvn.RData")
load("git_ALL_mvn230813.RData")

nilist = unlist(lapply(allgit, length))

pdat = data.frame( ID = rep(idvec, nilist),
                   index = unlist(lapply(nilist, function(x) 1:x)),
                   git = unlist(allgit), time = refdat$YearsSinceOnset )
thetam = apply(alltheta, 1, mean)
tmp = data.frame( ID = idvec, index = 0, git = thetam, time = 0)
pdat = rbind(tmp, pdat)
pdat = pdat[order(pdat$ID, pdat$index),]


#https://rpubs.com/droach/CPP526-codethrough

tmp = expand.grid(unique(pdat$ID), c(2,4,6,8,10))
tmp1 = data.frame(ID = tmp[,1], index = 999, git = NA, time = tmp[,2])
ppdat = rbind(pdat, tmp1)
ppdat = ppdat[order(ppdat$ID, ppdat$time),]

gap = 1
for(i in 1:nrow(ppdat)){
  if(ppdat$index[i] == 999){
    if(i == nrow(ppdat)){
      t = ppdat$time[i]
      dprev = 999
      if(ppdat$ID[i] == ppdat$ID[i-1]){
        tprev = ppdat$time[i-1];dprev = t - tprev
      }
      mgap = dprev
      if( mgap < gap){
        ppdat$git[i] = ppdat$git[i-1]
      }
    } else {
      t = ppdat$time[i]
      dprev = dnext = 999
      if(ppdat$ID[i] == ppdat$ID[i-1]){
        tprev = ppdat$time[i-1];dprev = t - tprev
      }
      if(ppdat$ID[i] == ppdat$ID[i+1]){
        tnext = ppdat$time[i+1]; dnext = tnext - t
      }
      mgap = min(c(dprev, dnext))
      if( mgap < gap){
        if(mgap == dprev) ppdat$git[i] = ppdat$git[i-1]
        if(mgap == dnext) ppdat$git[i] = ppdat$git[i+1]
      }
    }
  }
}
plotd = ppdat[ppdat$index %in% c(0,999),c(1,3,4)]
plotd$gitcat[plotd$git > 0 & plotd$git <= 0.2 ] = "0-.2"
plotd$gitcat[plotd$git > 0.2 & plotd$git <= 0.4 ] = ".2-.4"
plotd$gitcat[plotd$git > 0.4 & plotd$git <= 0.6 ] = ".4-.6"
plotd$gitcat[plotd$git > 0.6 & plotd$git <= 0.8 ] = ".6-.8"
plotd$gitcat[plotd$git > 0.8 & plotd$git <= 1 ] = ".8-1"


save(plotd, file  = "./BMC submission 2024 March/sankey_plotd.RData")

####################################################
load("./BMC submission 2024 March/sankey_plotd.RData")

library(ggalluvial)
library(ggplot2)
library(tidyr)
library(dplyr)
head(plotd)

plotd = plotd[plotd$time <= 8,]
plotd <- plotd %>% mutate(Probability = factor(gitcat, levels = c(".8-1", ".6-.8", ".4-.6", ".2-.4", "0-.2")))
colnames(plotd)[3] = "Time"

p.sankey <- ggplot(plotd, aes(x = Time, y = 1, stratum = Probability, alluvium = ID, fill = Probability, group = Probability)) +
 # scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  scale_fill_brewer(type = "qual", palette = "Set2")+
  #geom_text(stat = "stratum", size = 3) +
  #theme(legend.position = "none") +
  ggtitle("Change in P(Fast progressor | data)") + labs(y = "Number of Patients", x = "Years since onset") +
  theme_minimal()

pdf(height = 5, width = 7, file = "./BMC submission 2024 March/sankeyplot.pdf")
p.sankey
dev.off()



###########################################################################
setDT(pdat)
fst = pdat[, git[1], by = ID]
lst = pdat[, git[.N], by = ID]
hist(lst$V1 - fst$V1)
tmp = cbind(fst, lst$V1, lst$V1 - fst$V1)
View(tmp[abs(tmp$V3)>0.5,])
View(pdat[ID==278,])
###########################################################################


