
setwd("S://TrajectoryAnalysis//STAN")

### Prep
library(matrixcalc)
library(dplyr)
library(codetools)
#library(rstan)
library(splines)
library(data.table)
library(mvtnorm)
library(rpart)
#library(rpart.plot)
#library(randomForest)
#library(randomForestExplainer)
#library(glmnet)
library(coefplot)
library(xtable)


nit = 5000
nwup = 2000
seednum = 1
nc = 10 # number of chains

min_Ni = 2 # at least min_Ni obs for each person
subsetN = F # take a subset of patients

cut_year <- NA

source("dataprocessing3_0-10_scl70.R")
#load("git_ALL_mvn.RData")
load("git_ALL_mvn230813.RData")

load("S://TrajectoryAnalysis//lung function measures//GLI 2022 comparisons//demo_240419.Rdata")


load("candat_baseline_scl70.Rdata")

# baseline pft
baselinepft <- refdat %>% arrange(pt_id, date) %>% filter(!duplicated(pt_id))

# grouping
pdat <- refdat
pdat$pit <- unlist(allgit)
gdat <- pdat %>% arrange(desc(date)) %>% filter(!duplicated(pt_id))
gdat <- gdat %>% mutate(group = ifelse(pit > 0.5, "fast progressor" , "stable"))


data <- demo %>% merge(baselinepft, by = "pt_id") %>%

  mutate(onsetdate = pmin(rp_dt, nonrp_dt),

         ageonset_rp = as.numeric(rp_dt - DOB)/365.25,
         ageonset_nonrp = as.numeric(nonrp_dt - DOB)/365.25,
         ageonset = as.numeric(onsetdate - DOB)/365.25,

         diseasedur_rp = as.numeric(date - rp_dt)/365.25,
         diseasedur_nonrp = as.numeric(date - nonrp_dt)/365.25,
         diseasedur = as.numeric(date - onsetdate)/365.25,

         ageatbaseline = as.numeric(date - DOB)/365.25) %>%

  merge(gdat %>% select(pt_id, pit, group))

#data <- cbind(data, demodat %>% select(subtype, ageonset_rp, ageonset_nonrp, ageonset,
#                                       diseasedur_rp, diseasedur_nonrp, diseasedur))

mean_sd <- function(varname){
  vardat <- data %>% select(all_of(varname)); colnames(vardat) <- "vname"
  vtab <- vardat %>%
    summarise(g_mean = mean(vname, na.rm = T),
              g_sd = sd(vname, na.rm = T),
              N = sum(!is.na(vname)),
              perc_mis = round(sum(is.na(vname)) / length(vname )*100,2) )
  vtab <- vtab %>% mutate(vres = paste0(round(g_mean, 1), " (", round(g_sd, 1), "), ", perc_mis))
  return(res = vtab$vres)
}

mean_sd_noN <- function(varname){
  vardat <- data %>% select(all_of(varname)); colnames(vardat) <- "vname"
  vtab <- vardat %>%
    summarise(g_mean = mean(vname, na.rm = T),
              g_sd = sd(vname, na.rm = T),
              N = sum(!is.na(vname)))
  vtab <- vtab %>% mutate(vres = paste0(round(g_mean, 1), " (", round(g_sd, 1), ")"))
  return(res = vtab$vres)
}

mean_sd_group <- function(varname){
  vardat <- data %>% select(all_of(varname)); colnames(vardat) <- "vname"
  vtab <- vardat %>% group_by(data$group) %>%
    summarise(g_mean = mean(vname, na.rm = T),
              g_sd = sd(vname, na.rm = T),
              N = sum(!is.na(vname)),
              perc_mis = round(sum(is.na(vname)) / length(vname )*100,2) )
  vtab <- vtab %>% mutate(vres = paste0(round(g_mean, 1), " (", round(g_sd, 1), "), ", perc_mis))
  return(res = vtab$vres)
}


mean_sd_group_noN <- function(varname){
  vardat <- data %>% select(all_of(varname)); colnames(vardat) <- "vname"
  vtab <- vardat %>% group_by(data$group) %>%
    summarise(g_mean = mean(vname, na.rm = T),
              g_sd = sd(vname, na.rm = T),
              N = sum(!is.na(vname)) )
  vtab <- vtab %>% mutate(vres = paste0(round(g_mean, 1), " (", round(g_sd, 1), ")"))
  return(res = vtab$vres)
}

# Table 1: description of SSc cohort

t1demo <- rbind(mean_sd_noN("ageonset_nonrp"),
                mean_sd("ageonset_rp"),
                mean_sd("ageonset"),
                mean_sd("diseasedur_nonrp"),
                mean_sd("diseasedur_rp"),
                mean_sd("diseasedur"))

rownames(t1demo) <- c("ageonset_nonrp", "ageonset_rp", "ageonset",
                      "diseasedur_nonrp", "diseasedur_rp", "diseasedur")

t1pft <- rbind(mean_sd("ageatbaseline"),
               mean_sd_noN("fvcorig"),
               mean_sd_noN("dlcoorig"))

rownames(t1pft) <- c("ageatbaseline", "fvcorig", "dlcoorig")


# SSc Type (sine, limited, diffuse)
stab <- data %>%
  summarise(N = sum(!is.na(subtype)),
            N_sine = sum(subtype == "Sine", na.rm = T),
            perc_sine = N_sine/N * 100,
            N_limited = sum(subtype %in% c("Limited Type 1", "Limited Type 2"), na.rm = T),
            perc_limited = N_limited/N * 100,
            N_diffuse = sum(subtype == "Diffuse", na.rm = T),
            perc_diffuse = N_diffuse/N * 100)

t1type <- stab %>% mutate(sine = paste0(N_sine, " (", round(perc_sine, 1), "), N=", N),
                limited = paste0(N_limited, " (", round(perc_limited, 1), "), N=", N),
                diffuse = paste0(N_diffuse, " (", round(perc_diffuse, 1), ")")) %>%
  select(sine, limited, diffuse) %>% t

# sex
gtab <- data %>%
  summarise(N = sum(!is.na(male)),
            N_male = sum(male == 1, na.rm = T),
            perc_male = N_male/N * 100,
            N_female = sum(male == 0, na.rm = T),
            perc_female = N_female/N * 100)

t1sex <- gtab %>% mutate(male = paste0(N_male, " (", round(perc_male, 1), ")"),
                female = paste0(N_female, " (", round(perc_female, 1), ")")) %>%
  select(female, male) %>% t

# race ethnicity
gtab <- data %>%
  summarise(N = sum(!is.na(AArace)),
            N_race = sum(AArace == 1, na.rm = T),
            perc_race = N_race/N * 100)

t1race <- gtab %>% mutate(AArace = paste0(N_race, " (", round(perc_race, 1), ")")) %>%
  select(AArace) %>% t


# table by group
t2demo <- rbind(mean_sd_group_noN("ageonset_nonrp"),
                mean_sd_group("ageonset_rp"),
                mean_sd_group("ageonset"),
                mean_sd_group("diseasedur_nonrp"),
                mean_sd_group("diseasedur_rp"),
                mean_sd_group("diseasedur"))

rownames(t2demo) <- c("ageonset_nonrp", "ageonset_rp", "ageonset",
                      "diseasedur_nonrp", "diseasedur_rp", "diseasedur")

t2pft <- rbind(mean_sd_group("ageatbaseline"),
               mean_sd_group_noN("fvcorig"),
               mean_sd_group_noN("dlcoorig"))

rownames(t2pft) <- c("ageatbaseline", "fvcorig", "dlcoorig")

# SSc Type (sine, limited, diffuse)
stab <- data %>% group_by(data$group) %>%
  summarise(N = sum(!is.na(subtype)),
            N_sine = sum(subtype == "Sine", na.rm = T),
            perc_sine = N_sine/N * 100,
            N_limited = sum(subtype %in% c("Limited Type 1", "Limited Type 2"), na.rm = T),
            perc_limited = N_limited/N * 100,
            N_diffuse = sum(subtype == "Diffuse", na.rm = T),
            perc_diffuse = N_diffuse/N * 100)

t2type <- stab %>% mutate(sine = paste0(N_sine, " (", round(perc_sine, 1), "), N=", N),
                limited = paste0(N_limited, " (", round(perc_limited, 1), "), N=", N),
                diffuse = paste0(N_diffuse, " (", round(perc_diffuse, 1), ")")) %>%
  select(sine, limited, diffuse) %>% t

# sex
gtab <- data %>% group_by(data$group) %>%
  summarise(N = sum(!is.na(male)),
            N_male = sum(male == 1, na.rm = T),
            perc_male = N_male/N * 100,
            N_female = sum(male == 0, na.rm = T),
            perc_female = N_female/N * 100)

#t2sex <- gtab %>% mutate(male = paste0(N_male, " (", round(perc_male, 1), "), N=", N),
#                         female = paste0(N_female, " (", round(perc_female, 1), "), N=", N)) %>%
#  select(female, male) %>% t

t2sex <- gtab %>% mutate(male = paste0(N_male, " (", round(perc_male, 1), ")"),
                female = paste0(N_female, " (", round(perc_female, 1), ")")) %>%
  select(female, male) %>% t

# race ethnicity
gtab <- data %>% group_by(data$group) %>%
  summarise(N = sum(!is.na(AArace)),
            N_race = sum(AArace == 1, na.rm = T),
            perc_race = N_race/N * 100)

t2race <- gtab %>% mutate(AArace = paste0(N_race, " (", round(perc_race, 1), ")")) %>%
  select(AArace) %>% t

# AA (Euroimmun)
load("S://TrajectoryAnalysis//read in PMAP data//ebo_wide.Rdata")

AAdata <- merge(data, ebo_wide %>% select(-aca, -rnapol, -scl70), by = "pt_id", all.x = T)

summary_binvar <- function(dat, varname){
  vardat <- dat %>% select(all_of(varname)); colnames(vardat) <- "vname"
  vtab <- vardat %>%
    summarise(Npos = sum(vname == "1", na.rm = T),
              N = sum(!is.na(vname)),
              perc_pos = Npos/N * 100,
              perc_mis = round(sum(is.na(vname)) / length(vname )*100,2) )
  vtab <- vtab %>% mutate(vres = paste0(Npos, " (", round(perc_pos, 1), "), ", perc_mis))
  return(res = vtab$vres)
}

summary_binvar_group <- function(dat, varname){
  vardat <- dat %>% select(all_of(varname)); colnames(vardat) <- "vname"
  vtab <- vardat %>% group_by(AAdata$group) %>%
    summarise(Npos = sum(vname == "1", na.rm = T),
              N = sum(!is.na(vname)),
              perc_pos = Npos/N * 100,
              perc_mis = round(sum(is.na(vname)) / length(vname )*100,2) )
  vtab <- vtab %>% mutate(vres = paste0(Npos, " (", round(perc_pos, 1), "), ", perc_mis))
  return(res = vtab$vres)
}

AAname <- c("aca", "fib", "ku", "nor90", "PMScl", "rnapol", "ro52", "thto", "u1rnp")

AAres <- unlist(lapply(AAname, function(x)(summary_binvar(dat = AAdata, varname = x))))
AAmat <- matrix(AAres, nrow = 1, ncol = length(AAname)) %>% t
rownames(AAmat) <- AAname

AAres_group <- unlist(lapply(AAname, function(x)(summary_binvar_group(dat = AAdata, varname = x))))
AAmat_group <- matrix(AAres_group, nrow = 2, ncol = length(AAname)) %>% t
rownames(AAmat_group) <- AAname


tab1 <- cbind(rbind(t1demo, t1type, t1sex, t1race, t1pft, AAmat),
              rbind(t2demo, t2type, t2sex, t2race, t2pft, AAmat_group))

colnames(tab1) <- c("All", "Fast progressor group", "Stable group")

tab1 = tab1[c(11,12,9,1,14,15,     2:6,13,16:24),]

# All 289, Fast progressor 139, Stable 150

save(tab1, file = "tab1.RData")
xtable(tab1)
