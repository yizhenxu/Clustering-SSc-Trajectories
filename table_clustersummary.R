### table 1: clinical and demographic characteristics of clusters

### pull out posterior probability of being in mixture group 1
dn = dimnames(res)
lp = res[, , which(dn$parameters == "lp__")]
odch = order(apply(lp, 2, max), decreasing = T)

od = 1; idod = 1
maxidx = order(lp[, odch[od]],decreasing=T)[idod]

logpyz = matrix(res[maxidx, odch[od], grep("logpyz",dn$parameters )],ncol=2)
postz1 = rep(NA, nrow(logpyz))
for (n in 1:nrow(logpyz)) {
  logdenom = log(exp(logpyz[n,1]) + exp(logpyz[n,2]))
  postz1[n] = exp(logpyz[n,1] - logdenom)
}

gdat <- data.frame(prob1 = postz1) %>%
  mutate(group = ifelse(postz1 > 0.5, "group 1", "group 2"))

groupdat <- cbind(gdat, refdat %>% filter(!duplicated(pt_id)))
group1 <- groupdat %>% filter(group == "group 1")
group2 <- groupdat %>% filter(group == "group 2")

gtab <- groupdat %>% group_by(group) %>% summarise(
  N = length(group),
  male = mean(male),
  AArace = mean(AArace),
  ageonset = mean(ageonset),
  diffuse = mean(diffuse, na.rm = T),
  aca = mean(as.numeric(aca), na.rm = T),
  rnapol = mean(as.numeric(rnapol), na.rm = T),
  scl70 = mean(as.numeric(scl70), na.rm = T)) %>%
  as.data.frame()

gtab[, -c(1, 2, 5)] <- round(gtab[, -c(1, 2, 5)]*100, 1)

gtab


#### read in data for ILD - pulmonary fibrosis ####

# packages
Packages <- c("DBI", "RODBC")
lapply(Packages, library, character.only = T)

# dbi connection
connectionString <- "Driver={SQL Server};server=ESMPMDBPR4.WIN.AD.JHU.EDU;database=Scleroderma_Projection;trusted_connection=Yes;"
odbccon <- odbcDriverConnect(connection = connectionString)

# EPIC MRN
sql <- "SELECT [cohort_id]
      ,[pt_id]
      ,[visit_dt]
      ,[visit_pulmfibrosis]
      ,[epic_ildcriteria]
	  FROM [Scleroderma_PHI_Projection].[dbo].[visit_assessment]"

ilddat <- sqlQuery(odbccon, sql)

ildid <- ilddat[ilddat$visit_pulmfibrosis == 1, "pt_id"] %>% unique # 1663

ilddat$cohort_id %>% unique %>% length # 4197

head(groupdat)

g1nILD <- groupdat[groupdat$group == "group 1", "pt_id"] %in% ildid %>% sum
g2nILD <- groupdat[groupdat$group == "group 2", "pt_id"] %in% ildid %>% sum


data.frame(gtab, nILD = c(g1nILD, g2nILD)) %>% mutate(pILD = round(nILD/N, 2))



