
##################################################

# 'ncores' is the number of cores to be used
# for parallel computing in Windows with 64-bit R:
ncores <- 10

# 'nsims' is the number of simulated data sets:
nsims <- 500

## Importing source code & loading packages:
library(MCMCpack)
library("parallel")

source("./source_code/MainFunctions.R")
source("./source_code/SupportPrograms.R")
source("./source_code/GT_sim_fns.R")

# Making clusters for parallel computing:
cl <- makeCluster(ncores)
clusterExport(cl, "hier.alg.data")
clusterExport(cl, "array.2dim.data")
clusterExport(cl, "PostUnknownAssayAcr")
clusterExport(cl, "poolMembTracker")
clusterExport(cl, "EmUnknownAssayAcr")
clusterExport(cl, "EmultGibbs2")
clusterExport(cl, "PostKnownAssayAcr")
clusterExport(cl, "EmKnownAssayAcr")
clusterExport(cl, "EmultGibbs")

##################################################

## The code below reproduces the simulation results 
## depicted in Table D.1 of the Supporting Information.
## The accuracy parameter 'delta' is assumed to be known.

# Setting seed:
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(123)
clusterSetRNGStream(cl = cl, iseed = 123)

# Simulation parameters:
N <- 5000         # Sample size
Se <- c(.95,.95)  # True sensitivities for diseases 1 & 2
Sp <- c(.99,.99)  # True specificities for diseases 1 & 2

## Specify the prior hyperparameters as shown below. 
## Flat priors have been used for all parameters.
p.pr <- rep(1,4)
Se1.pr <- c(1,1)
Se2.pr <- c(1,1)
Sp1.pr <- c(1,1)
Sp2.pr <- c(1,1)

G <- 12000               # number of Gibbs iterates
burn <- 2000             # burn-in period
pick <- seq(1,10000,5)   # keeping every 5th

## Convergence tolerance for the EM algorithm
epsilon <- 0.001

## Maximum number of iterations for the EM algorithm. 
emmaxit <- 200

## Choose an initial value of the parameter p. The initial value 
## can be specified from historical, pilot study, or any contextual 
## estimate. We use p0 <- c(.92,.05,.02,.01) throughout the paper.
p0 <- c(.92,.05,.02,.01)

options( scipen = 999 )
Table.D1 <- matrix(NA, nrow = 25, ncol = 19)

for(Configuration in c("I", "II")){

if(Configuration == "I"){
  ## Configuration I
  p <- c(.95,.02,.02,.01)     # True value of p
  MPT <- 5
  H2 <- c(5,1)
  H3 <- c(9,3,1)
  H4 <- c(18,6,3,1)
  AT <- c(11,11,1)  # For AT, c(11,11,1), not c(11,11)
                    # because tests are done in 3 stages:
                    # row pool, col pool, & individual testing
  i1 <- 0
  i2 <- 0
}

if(Configuration == "II"){
  ## Configuration II
  p <- c(.990,.004,.004,.002) # True value of p
  MPT <- 11
  H2 <- c(11,1)
  H3 <- c(25,5,1)
  H4 <- c(48,12,4,1)
  AT <- c(29,29,1)
  i1 <- 13
  i2 <- 10
}

#========= Estimation with MPT protocol =========#

protocol <- MPT
list.zmat <- lapply( rep(N, nsims), hier.alg.data, p=p,
                     design=protocol, Se=Se, Sp=Sp )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, N=N, S=length(protocol), p.pr=p.pr,
           postGit=G, method="Bayes", accuracy="known")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, N=N, S=length(protocol), p.pr=p.pr,
           emGit=G, emburn=burn, method="MAP", accuracy="known")


p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
T <- rep(-9, nsims)

for(s in 1:nsims){
## POSTERIOR SAMPLING 
  res1 <- res.bayes[[s]]
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D1[i1+1:4, 1] <- round(colMeans(p.Mean),3)   # Est
Table.D1[i1+1:4, 2] <- round(apply(p.Mean,2,sd),4) # SD
Table.D1[i1+1:4, 3] <- round(colMeans(p.SE),4)     # SE

## MAP Estimation:
Table.D1[i1+6:9, 1] <- round(colMeans(p.MAP),3)    # Est
Table.D1[i1+6:9, 2] <- round(apply(p.MAP,2,sd),4)  # SD
Table.D1[i1+6:9, 3] <- round(colMeans(p.SE),4)     # SE

## Avg. number of tests:
Table.D1[i1+11, 2] <- round( mean(T), 1 )    
Table.D1[i1+12, 2] <- round( sd(T), 1 )


#========= Estimation with H2 protocol =========#

protocol <- H2
list.zmat <- lapply( rep(N, nsims), hier.alg.data, p=p,
                     design=protocol, Se=Se, Sp=Sp )
		   
## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, N=N, S=length(protocol), p.pr=p.pr,
           postGit=G, method="Bayes", accuracy="known")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, N=N, S=length(protocol), p.pr=p.pr,
           emGit=G, emburn=burn, method="MAP", accuracy="known")

p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
T <- rep(-9, nsims)

for(s in 1:nsims){
## POSTERIOR SAMPLING 
  res1 <- res.bayes[[s]]
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D1[i1+1:4, 5] <- round(colMeans(p.Mean),3)   # Est
Table.D1[i1+1:4, 6] <- round(apply(p.Mean,2,sd),4) # SD
Table.D1[i1+1:4, 7] <- round(colMeans(p.SE),4)     # SE

## MAP Estimation:
Table.D1[i1+6:9, 5] <- round(colMeans(p.MAP),3)    # Est
Table.D1[i1+6:9, 6] <- round(apply(p.MAP,2,sd),4)  # SD
Table.D1[i1+6:9, 7] <- round(colMeans(p.SE),4)     # SE

## Avg. number of tests:
Table.D1[i1+11, 6] <- round( mean(T), 1 )    
Table.D1[i1+12, 6] <- round( sd(T), 1 )

#========= Estimation with H3 protocol =========#

protocol <- H3
list.zmat <- lapply( rep(N, nsims), hier.alg.data, p=p,
                     design=protocol, Se=Se, Sp=Sp )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, N=N, S=length(protocol), p.pr=p.pr,
           postGit=G, method="Bayes", accuracy="known")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, N=N, S=length(protocol), p.pr=p.pr,
           emGit=G, emburn=burn, method="MAP", accuracy="known")

p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
T <- rep(-9, nsims)

for(s in 1:nsims){
## POSTERIOR SAMPLING 
  res1 <- res.bayes[[s]]
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D1[i1+1:4, 9] <- round(colMeans(p.Mean),3)   # Est
Table.D1[i1+1:4, 10] <- round(apply(p.Mean,2,sd),4) # SD
Table.D1[i1+1:4, 11] <- round(colMeans(p.SE),4)     # SE

## MAP Estimation:
Table.D1[i1+6:9, 9] <- round(colMeans(p.MAP),3)    # Est
Table.D1[i1+6:9, 10] <- round(apply(p.MAP,2,sd),4)  # SD
Table.D1[i1+6:9, 11] <- round(colMeans(p.SE),4)     # SE

## Avg. number of tests:
Table.D1[i1+11, 10] <- round( mean(T), 1 )    
Table.D1[i1+12, 10] <- round( sd(T), 1 )

#========= Estimation with H4 protocol =========#

protocol <- H4
list.zmat <- lapply( rep(N, nsims), hier.alg.data, p=p,
                     design=protocol, Se=Se, Sp=Sp )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, N=N, S=length(protocol), p.pr=p.pr,
           postGit=G, method="Bayes", accuracy="known")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, N=N, S=length(protocol), p.pr=p.pr,
           emGit=G, emburn=burn, method="MAP", accuracy="known")

p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
T <- rep(-9, nsims)

for(s in 1:nsims){
## POSTERIOR SAMPLING 
  res1 <- res.bayes[[s]]
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D1[i1+1:4, 13] <- round(colMeans(p.Mean),3)   # Est
Table.D1[i1+1:4, 14] <- round(apply(p.Mean,2,sd),4) # SD
Table.D1[i1+1:4, 15] <- round(colMeans(p.SE),4)     # SE

## MAP Estimation:
Table.D1[i1+6:9, 13] <- round(colMeans(p.MAP),3)    # Est
Table.D1[i1+6:9, 14] <- round(apply(p.MAP,2,sd),4)  # SD
Table.D1[i1+6:9, 15] <- round(colMeans(p.SE),4)     # SE

## Avg. number of tests:
Table.D1[i1+11, 14] <- round( mean(T), 1 )    
Table.D1[i1+12, 14] <- round( sd(T), 1 )

#========= Estimation with AT protocol =========#

protocol <- AT
list.zmat <- lapply( rep(N, nsims), array.2dim.data, p=p,
                     design=protocol, Se=Se, Sp=Sp )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, N=N, S=length(protocol), p.pr=p.pr,
           postGit=G, method="Bayes", accuracy="known")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, N=N, S=length(protocol), p.pr=p.pr,
           emGit=G, emburn=burn, method="MAP", accuracy="known")

p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
T <- rep(-9, nsims)

for(s in 1:nsims){
## POSTERIOR SAMPLING 
  res1 <- res.bayes[[s]]
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D1[i1+1:4, 17] <- round(colMeans(p.Mean),3)   # Est
Table.D1[i1+1:4, 18] <- round(apply(p.Mean,2,sd),4) # SD
Table.D1[i1+1:4, 19] <- round(colMeans(p.SE),4)     # SE

## MAP Estimation:
Table.D1[i1+6:9, 17] <- round(colMeans(p.MAP),3)    # Est
Table.D1[i1+6:9, 18] <- round(apply(p.MAP,2,sd),4)  # SD
Table.D1[i1+6:9, 19] <- round(colMeans(p.SE),4)     # SE

## Avg. number of tests:
Table.D1[i1+11, 18] <- round( mean(T), 1 )    
Table.D1[i1+12, 18] <- round( sd(T), 1 )

}

stopCluster( cl )

## Preparing & saving Table D.1
dfr1 <- data.frame( cbind(
     c("","","","","","Config-I","","","","","","", 
	   "",  "","","","","","Config-II","","","","","",""),
     c("", "", "", "p00=0.95", "p10=0.02", "p01=0.02", "p11=0.01", "", "", "", "", "", "", 
	   "", "", "", "p00=0.990", "p10=0.004", "p01=0.004", "p11=0.002", "", "", "", "", ""),  
     c("Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T", "", 
	   "Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T")))

dfr1 <- data.frame( dfr1, Table.D1 )
dfr1 <- sapply(dfr1, as.character)
dfr1[is.na(dfr1)] <- " "
Table_D1 <- as.data.frame(dfr1)
colnames( Table_D1 ) <- c("", "","", "Est", "SD", "SE", "", "Est", "SD", "SE", 
                        "", "Est", "SD", "SE", "", "Est", "SD", "SE", "", "Est", "SD", "SE")

saveRDS(Table_D1, "./simulation/intermediate_results/Table_D1.rds")

