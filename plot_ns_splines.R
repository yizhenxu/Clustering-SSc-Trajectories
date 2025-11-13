
setwd("S:/TrajectoryAnalysis/STAN")

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

source("dataprocessing3_0-10_scl70.R")

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

tx = seq(0.01,9.99,0.01)
tmp = ns(tx, knots = 5, Boundary.knots = c(0, 10))


pdf("./BMC submission 2024 March/Supp_splines.pdf", width = 6, height = 3)
par(mfrow = c(1,2),mar = c(4,4,1,1))
plot(tx, tmp[,1], type = "l", xlab = "Years since onset", ylab = "Spline 1")
plot(tx, tmp[,2], type = "l", xlab = "Years since onset", ylab = "Spline 2")
dev.off()

