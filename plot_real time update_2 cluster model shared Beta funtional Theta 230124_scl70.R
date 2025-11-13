
setwd("S://TrajectoryAnalysis//STAN")

# notes for Yizhen: script with new function gij_mvn

# load results: beta1, beta2, betaS, G1, G2, R
source("plot_table_preparation.R") # related to task 2

# load functions: git_fun
source("Functions.R")

### Prep
library(matrixcalc)
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
### wrap up a new function as function of pt_id, for all the 289 ppl

library(parallel)
detectCores() # 6 cores

wrapupj_showargs = function(j, Yalli, xYi, xSi, Zi, beta1, beta2, betaS, G1, G2,R,alpha){


  beta1j = t(matrix(beta1[j,], nrow = 2))
  beta2j = t(matrix(beta2[j,], nrow = 2))
  betaSj = t(matrix(betaS[j,], nrow = 2))
  G1j = matrix(G1[j, ], nrow = 2)
  G2j = matrix(G2[j, ], nrow = 2)
  Rj = matrix(R[j, ], nrow = 2)
  alphaj = alpha[j, ]

  res = gij_mvn(Yalli, xYi, xSi, beta1j, beta2j, betaSj, G1j, G2j, Rj, Zi, alphaj)

  return(t(res))

}

wrapup_all = function(id, refdat, xY, xS, Z, beta1, beta2, betaS, G1, G2,R,alpha){

  pidind <- which(refdat$pt_id == id)
  Yalli = refdat[pidind, list(fvc,dlco)]
  xYi <- xY[pidind, ]
  xSi <- xS[pidind, ]
  Zi <- Z[pidind, ]

  git_mat = simplify2array(lapply(1:npost, function(j) wrapupj_showargs(j, Yalli, xYi, xSi, Zi,beta1, beta2, betaS, G1, G2, R, alpha)))
  git_mat = git_mat[1,,]
  # posterior mean
  git = as.numeric(apply(git_mat, 1, mean))

  return(git)
}

npost = dim(lp)[1]
person_i =  function(id) wrapup_all(id, refdat, xY, xS, Z, beta1, beta2, betaS, G1, G2, R, alpha)
#person_i(20)

#idvec = unique(refdat$pt_id)[1:3] # test for the first 3 ppl
idvec = unique(refdat$pt_id)


start = Sys.time()
allgit = mclapply(idvec, person_i)
end = Sys.time()
end - start # for 6 ppl, took Time difference of 3.171589 mins


#save(allgit, file  = "git_ALL_mvn.RData")

save(allgit, file  = "git_ALL_mvn230813.RData") # 3h


###########################################################################
# calculating theta = P(cluster 1)

#load("git_ALL_mvn.RData")
load("git_ALL_mvn230813.RData")

Ztmp1 <- refdat %>% filter(!duplicated(pt_id)) %>% select(all_of(factorname))
Ztmp2 = cbind(1, Ztmp1)

# matrix of N x npost  (289 x 3000)
library(boot) # inv.logit
alltheta = inv.logit(as.matrix(Ztmp2) %*% t(alpha))
#hist(alltheta)

###########################################################################

# number of patients in each cluster
if(0){
  pvecall <- unlist(lapply(allgit, function(x)tail(x, n = 1)))
  sum(pvecall > 0.5)
}

nilist = unlist(lapply(allgit, length))

pdat = data.frame( ID = rep(idvec, nilist),
                   index = unlist(lapply(nilist, function(x) 1:x)),
                   git = unlist(allgit), time = refdat$YearsSinceOnset )
thetam = apply(alltheta, 1, mean)
tmp = data.frame( ID = idvec, index = 0, git = thetam, time = 0)
pdat = rbind(tmp, pdat)
pdat = pdat[order(pdat$ID, pdat$index),]

if(0){
  gdat <- pdat %>% arrange(ID, time) %>% group_by(ID) %>% summarise(lastp = last(git)) %>%
    mutate(group = ifelse(lastp < 0.5, "progressor", "stable"))

  save(gdat, file = "ID_by_group_scl70.Rdata")
}

#png("./res/plot_ptheta_spaghetti.png", width = 700, height = 400)

ptab <- pdat[!duplicated(pdat$ID), c("ID", "git")]
pdat <- merge(pdat, ptab, by = "ID", all.x = T)
colnames(pdat)[5] <- "Probability"

library(RColorBrewer)
mypalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_color_gradientn(colours = mypalette(100), limits = c(0, 1))

pdf(file = "./BMC submission 2024 March/clusterprob_all.pdf", width = 8, height = 6)

ggplot(pdat, aes(x=time, y=git.x, group = ID, color = Probability)) +
  geom_line() + sc +
  #guides(colour="none") +
  xlab("Years Since Onset") +
  ylab("Probability in Cluster 1")#+ aes(colour = factor(ID))

dev.off()

if(0){
  setDT(pdat)
  fst = pdat[, git[1], by = ID]
  lst = pdat[, git[.N], by = ID]
  hist(lst$V1 - fst$V1)
  tmp = cbind(fst, lst$V1, lst$V1 - fst$V1)
  View(tmp[abs(tmp$V3)>0.5,])
  View(pdat[ID==278,])
}


refdat = as.data.table(refdat)
idvec = unique(refdat$pt_id)

library(tidyr)

id = 487#3673
dn = dimnames(res)
lp = res[, , which(dn$parameters == "lp__")]
npost = dim(lp)[1]
# extreme IDs: 278, 487, 516, 614, 670, 795, 1109, 1250, 2514, 3646, 3673

odch = order(apply(lp, 2, mean), decreasing = T)
alldrawslst = lapply(sample(size = 1000, 1:npost), function(j) eachperson_NP_xS_pred(id, odch[1], j))

alldraws = simplify2array(alldrawslst)  # array of Ti x 4 x npost

# pdat
pdatind <- pdat %>% filter(ID == id)

nT = nrow(refdat[pt_id == id,])

pd = data.frame(time = rep(refdat[pt_id == id,YearsSinceOnset], 4),
                fvc = rep(refdat[pt_id == id, fvc], 4),
                dlco = rep(refdat[pt_id == id, dlco], 4),
                fvcorig = rep(refdat[pt_id == id, fvcorig], 4),
                dlcoorig = rep(refdat[pt_id == id, dlcoorig], 4),

                m = c(apply(alldraws, c(1,2), mean)),
                l = c(apply(alldraws, c(1,2), lfun)),
                u = c(apply(alldraws, c(1,2), ufun)),

                type = rep(c("FVC1", "DLCO1", "FVC2", "DLCO2"), each=nT))

obsdat <- pd %>% select(time, fvc, dlco) %>% gather(key = measure, value = value, fvc:dlco)


maxtime <- max(pdatind$time)

library(ggplot2)

p1 <- ggplot()+
  geom_line(data = pd %>% filter(type %in% c("FVC1", "FVC2")),
            aes(x = time, y = m, group = type, color = type))+
  geom_ribbon(data = pd %>% filter(type %in% c("FVC1", "FVC2")),
              aes(x = time, y = m, group = type,
                  fill = as.factor(type), ymin = l, ymax=u), alpha=0.3, show.legend = F) +

  geom_line(data = obsdat %>% filter(measure == "fvc"), aes(x = time, y = value)) +
  geom_point(data = obsdat %>% filter(measure == "fvc"), aes(x = time, y = value)) +
  xlim(0, maxtime) + theme(legend.position = "right") + labs(y = "pFVC", x = "", color = "") +
  scale_color_discrete(labels = c("Fast progressor", "Stable"))

p2 <- ggplot()+
  geom_line(data = pd %>% filter(type %in% c("DLCO1", "DLCO2")),
            aes(x = time, y = m, group = type, color = type))+
  geom_ribbon(data = pd %>% filter(type %in% c("DLCO1", "DLCO2")),
              aes(x = time, y = m, group = type,
                  fill = as.factor(type), ymin = l, ymax=u), alpha=0.3) +

  geom_line(data = obsdat %>% filter(measure == "dlco"), aes(x = time, y = value)) +
  geom_point(data = obsdat %>% filter(measure == "dlco"), aes(x = time, y = value)) +
  xlim(0, maxtime) + theme(legend.position = "none") + labs(y = "pDLCO", x = "")

p3 <- ggplot() + geom_point(data = pdatind, aes(x = time, y = git.x)) +
  geom_line(data = pdatind, aes(x = time, y = git.x)) + labs(y = "P(Fast progressor | data)",x="Years since onset") +
  ylim(0, 1) + xlim(0, maxtime)


pdf(height = 7, width = 5, file = paste0("./BMC submission 2024 March/Figure4_Predtraj_pid", id,"_240709.pdf"))
ggpubr::ggarrange(p1, p2, p3, nrow = 3, common.legend = T)
dev.off()

# time -> years since onset
# FVCt, DLCTt -> Standardized ppFVC, Standardized ppDLCO, describe pp in caption
# check clusters 1,2 labels -> Fast progressor, Stable
# P(cluster 1|data) -> P(Fast progressor | data)

# 
# ### for loop for better life
# forloopids = c(516,278,487,29,39) #3673
# for(id in forloopids){
#   dn = dimnames(res)
#   lp = res[, , which(dn$parameters == "lp__")]
#   npost = dim(lp)[1]
#   # extreme IDs: 278, 487, 516, 614, 670, 795, 1109, 1250, 2514, 3646, 3673
#   
#   odch = order(apply(lp, 2, mean), decreasing = T)
#   alldrawslst = lapply(sample(size = 1000, 1:npost), function(j) eachperson_NP_xS_pred(id, odch[1], j))
#   
#   alldraws = simplify2array(alldrawslst)  # array of Ti x 4 x npost
#   
#   # pdat
#   pdatind <- pdat %>% filter(ID == id)
#   
#   nT = nrow(refdat[pt_id == id,])
#   
#   pd = data.frame(time = rep(refdat[pt_id == id,YearsSinceOnset], 4),
#                   fvc = rep(refdat[pt_id == id, fvc], 4),
#                   dlco = rep(refdat[pt_id == id, dlco], 4),
#                   fvcorig = rep(refdat[pt_id == id, fvcorig], 4),
#                   dlcoorig = rep(refdat[pt_id == id, dlcoorig], 4),
#                   
#                   m = c(apply(alldraws, c(1,2), mean)),
#                   l = c(apply(alldraws, c(1,2), lfun)),
#                   u = c(apply(alldraws, c(1,2), ufun)),
#                   
#                   type = rep(c("FVC1", "DLCO1", "FVC2", "DLCO2"), each=nT))
#   
#   obsdat <- pd %>% select(time, fvc, dlco) %>% gather(key = measure, value = value, fvc:dlco)
#   
#   
#   maxtime <- max(pdatind$time)
#   
#   p1 <- ggplot()+
#     geom_line(data = pd %>% filter(type %in% c("FVC1", "FVC2")),
#               aes(x = time, y = m, group = type, color = type))+
#     geom_ribbon(data = pd %>% filter(type %in% c("FVC1", "FVC2")),
#                 aes(x = time, y = m, group = type,
#                     fill = as.factor(type), ymin = l, ymax=u), alpha=0.3, show.legend = F) +
#     
#     geom_line(data = obsdat %>% filter(measure == "fvc"), aes(x = time, y = value)) +
#     geom_point(data = obsdat %>% filter(measure == "fvc"), aes(x = time, y = value)) +
#     xlim(0, maxtime) + theme(legend.position = "right") + labs(y = "FVCt", x = "", color = "") +
#     scale_color_discrete(labels = c("Cluster 1", "Cluster 2"))
#   
#   p2 <- ggplot()+
#     geom_line(data = pd %>% filter(type %in% c("DLCO1", "DLCO2")),
#               aes(x = time, y = m, group = type, color = type))+
#     geom_ribbon(data = pd %>% filter(type %in% c("DLCO1", "DLCO2")),
#                 aes(x = time, y = m, group = type,
#                     fill = as.factor(type), ymin = l, ymax=u), alpha=0.3) +
#     
#     geom_line(data = obsdat %>% filter(measure == "dlco"), aes(x = time, y = value)) +
#     geom_point(data = obsdat %>% filter(measure == "dlco"), aes(x = time, y = value)) +
#     xlim(0, maxtime) + theme(legend.position = "none") + labs(y = "DLCOt", x = "")
#   
#   p3 <- ggplot() + geom_point(data = pdatind, aes(x = time, y = git.x)) +
#     geom_line(data = pdatind, aes(x = time, y = git.x)) + labs(y = "p(cluster1|data)") +
#     ylim(0, 1) + xlim(0, maxtime)
#   
#   
#   pdf(height = 7, width = 5, file = paste0("./BMC submission 2024 March/Predtraj_pid", id,"_240709.pdf"))
#   ggpubr::ggarrange(p1, p2, p3, nrow = 3, common.legend = T)
#   dev.off()
#   
# }
