# Functional theta
# two clusters, shared R, year 0 ~ 10, scl70 positive


#install.packages("rstan",lib = "C:/Program Files/RPackages")
#.libPaths()

setwd("S:/TrajectoryAnalysis/STAN")

setwd("/home/idies/workspace/Storage/jkim478/persistent/rstan/")

setwd("/home/idies/workspace/Storage/yxu143/persistent/Scleroderma clustering/")

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

stanfile <- "stan_2cluster_sharedR_nR_b0_sharedBeta_functionalTheta.R" # stan parvec datalist location
fn <- "mod2_functionalTheta.stan" # stan file name to save as
fitname <- "mod2fit_sharedBeta_230813_scl70_functionalTheta.RData" # fit saved into this
resname <- "mod2res_sharedBeta_230813_scl70_functionalTheta.RData" # res saved into this

###################################################

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
Zmat <- refdat %>% filter(!duplicated(pt_id)) %>% select(all_of(factorname))

Z = cbind(1, Zmat)
colnames(Z)[1] <- "Intercept"


### number of mixture components L
### mixture membership is by person, indexed by i
### generated quantitiy 'logpyz': Nsubs x L matrix. The (i,l)th element being
###   the posterior probability of person i being a member of the lth group

source(stanfile) # define stan code, datalist, parvec

# Delete existing mod.stan
if (file.exists(fn))  file.remove(fn)

# Write a new one
writeLines(scode, fn)
rstan:::rstudio_stanc(fn)


options(mc.cores = nc) #parallel::detectCores()
fit = stan(file = fn, data = datalist,
           pars = parvec,
           include = T, chains = nc, iter = nit, warmup = nwup, thin = 1, seed = seednum)

save(fit, file  = fitname)

res = extract(fit, permuted = F, inc_warmup = F)
save(res, file = resname)

# traceplot - log posterior
pdf("plotlp.pdf")
traceplot(fit, pars = "lp__", inc_warmup = F)
dev.off()

# beta - group 1
pdf("plotbeta1.pdf")
traceplot(fit, pars = "beta_1", inc_warmup = F)
dev.off()

# beta - group 2
pdf("plotbeta2.pdf")
traceplot(fit, pars = "beta_2", inc_warmup = F)
dev.off()

# mixing probabilities
pdf("plotalpha.pdf")
traceplot(fit, pars = "alpha", inc_warmup = F)
dev.off()
