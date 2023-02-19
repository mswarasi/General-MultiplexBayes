
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

## This code reproduces simulation results depicted 
## in Tables D.9-D.11 of the Supporting Information.

# Setting seed:
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(123)
clusterSetRNGStream(cl = cl, iseed = 123)

# Sample size:
N <- 1000         

## Configuration I:
p <- c(.95,.02,.02,.01)     # True value of p
MPT <- 5
H2 <- c(5, 1)
H3 <- c(9, 3, 1)
H4 <- c(18, 6, 3, 1)
AT <- c(11, 11, 1) 

## Specify the prior hyperparameters as shown below. 
## Flat priors have been used for all parameters.
p.pr <- rep(1, 4)
Se1.pr <- c(1, 1)
Se2.pr <- c(1, 1)
Sp1.pr <- c(1, 1)
Sp2.pr <- c(1, 1)

G <- 12000               # number of Gibbs iterates
burn <- 2000             # burn-in period
pick <- seq(1,10000,5)   # keeping every 5th

## Convergence tolerance for the EM algorithm
epsilon <- 0.001

## Maximum number of iterations for the EM algorithm. 
emmaxit <- 200

## Remarks about the initial values:
## --------------------------------
## The proposed estimation techniques require an initial value of 
## p and delta. These values can be specified from historical, pilot
## study, or any contextual estimate. We use p0=c(.92,.05,.02,.01) and
## delta0=c(.96,.96,.98,.98) throughout. Note that p0=c(.92,.05,.02,.01) 
## is comparable to the vector of coinfection probabilities in our Iowa 
## CT/NG data and delta0=c(.96,.96,.98,.98) is consistent with the 
## accuracy of the Aptima Combo 2 Assay used for CT/NG testing.
p0 <- c(.92, .05, .02, .01)
delta0 <- c(.96, .96, .98, .98)

options( scipen = 999 )


##################################################

## Table D.9 simulation assuming KNOWN accuracy:

Table.D9 <- matrix(NA, nrow = 25, ncol = 19)

for(Accuracy in c("high", "low")){

  if(Accuracy == "high"){
    Se <- c(.95, .95) 
    Sp <- c(.99, .99)  
    i1 <- 0
    i2 <- 0
  }

  if(Accuracy == "low"){
    Se <- c(.85, .85) 
    Sp <- c(.90, .90) 
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
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn), ][pick, ])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn), ][pick, ], 2, sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D9[i1+1:4, 1] <- round(colMeans(p.Mean), 3)   # Est
Table.D9[i1+1:4, 2] <- round(apply(p.Mean,2,sd), 4) # SD
Table.D9[i1+1:4, 3] <- round(colMeans(p.SE), 4)     # SE

## MAP Estimation:
Table.D9[i1+6:9, 1] <- round(colMeans(p.MAP), 3)    # Est
Table.D9[i1+6:9, 2] <- round(apply(p.MAP, 2, sd), 4)  # SD
Table.D9[i1+6:9, 3] <- round(colMeans(p.SE), 4)     # SE

## Avg. number of tests:
Table.D9[i1+11, 2] <- round( mean(T), 1 )    
Table.D9[i1+12, 2] <- round( sd(T), 1 )


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
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn), ][pick, ])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn), ][pick, ], 2, sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D9[i1+1:4, 5] <- round(colMeans(p.Mean), 3)     # Est
Table.D9[i1+1:4, 6] <- round(apply(p.Mean, 2, sd), 4) # SD
Table.D9[i1+1:4, 7] <- round(colMeans(p.SE), 4)       # SE

## MAP Estimation:
Table.D9[i1+6:9, 5] <- round(colMeans(p.MAP),  3)     # Est
Table.D9[i1+6:9, 6] <- round(apply(p.MAP, 2, sd), 4)  # SD
Table.D9[i1+6:9, 7] <- round(colMeans(p.SE), 4)       # SE

## Avg. number of tests:
Table.D9[i1+11, 6] <- round( mean(T), 1 )    
Table.D9[i1+12, 6] <- round( sd(T), 1 )

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
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn), ][pick, ])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn), ][pick, ], 2, sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D9[i1+1:4, 9] <- round(colMeans(p.Mean), 3)   # Est
Table.D9[i1+1:4, 10] <- round(apply(p.Mean, 2, sd), 4) # SD
Table.D9[i1+1:4, 11] <- round(colMeans(p.SE), 4)     # SE

## MAP Estimation:
Table.D9[i1+6:9, 9] <- round(colMeans(p.MAP), 3)    # Est
Table.D9[i1+6:9, 10] <- round(apply(p.MAP, 2, sd), 4)  # SD
Table.D9[i1+6:9, 11] <- round(colMeans(p.SE), 4)     # SE

## Avg. number of tests:
Table.D9[i1+11, 10] <- round( mean(T), 1 )    
Table.D9[i1+12, 10] <- round( sd(T), 1 )

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
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn), ][pick, ])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn), ][pick, ], 2, sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D9[i1+1:4, 13] <- round(colMeans(p.Mean), 3)     # Est
Table.D9[i1+1:4, 14] <- round(apply(p.Mean, 2, sd), 4) # SD
Table.D9[i1+1:4, 15] <- round(colMeans(p.SE), 4)       # SE

## MAP Estimation:
Table.D9[i1+6:9, 13] <- round(colMeans(p.MAP),   3)    # Est
Table.D9[i1+6:9, 14] <- round(apply(p.MAP, 2, sd), 4)  # SD
Table.D9[i1+6:9, 15] <- round(colMeans(p.SE), 4)       # SE

## Avg. number of tests:
Table.D9[i1+11, 14] <- round( mean(T), 1 )    
Table.D9[i1+12, 14] <- round( sd(T), 1 )

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
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn), ][pick, ])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn), ][pick, ], 2, sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D9[i1+1:4, 17] <- round(colMeans(p.Mean), 3)   # Est
Table.D9[i1+1:4, 18] <- round(apply(p.Mean, 2, sd), 4) # SD
Table.D9[i1+1:4, 19] <- round(colMeans(p.SE), 4)     # SE

## MAP Estimation:
Table.D9[i1+6:9, 17] <- round(colMeans(p.MAP), 3)      # Est
Table.D9[i1+6:9, 18] <- round(apply(p.MAP, 2, sd), 4)  # SD
Table.D9[i1+6:9, 19] <- round(colMeans(p.SE), 4)       # SE

## Avg. number of tests:
Table.D9[i1+11, 18] <- round( mean(T), 1 )    
Table.D9[i1+12, 18] <- round( sd(T), 1 )

}

## Preparing & saving Table D.9
dfr1 <- data.frame( cbind(
     c("","","","","","High accuracy","","","","","","", 
	   "",  "","","","","","Low accuracy","","","","","",""),
     c("", "", "", "p00=0.95", "p10=0.02", "p01=0.02", "p11=0.01", "", "", "", "", "", "", 
	   "", "", "", "p00=0.95", "p10=0.02", "p01=0.02", "p11=0.01", "", "", "", "", ""),  
     c("Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T", "", 
	   "Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T")))

dfr1 <- data.frame( dfr1, Table.D9 )
dfr1 <- sapply(dfr1, as.character)
dfr1[is.na(dfr1)] <- " "
Table_D9 <- as.data.frame(dfr1)
colnames( Table_D9 ) <- c("", "","", "Est", "SD", "SE", "", "Est", "SD", "SE", 
                        "", "Est", "SD", "SE", "", "Est", "SD", "SE", "", "Est", "SD", "SE")

saveRDS(Table_D9, "./simulation/intermediate_results/Table_D9.rds")


##################################################

## Table D.10-D.11 simulation assuming UNKNOWN accuracy:

Table.D10 <- matrix(NA, nrow = 25, ncol = 15)
Table.D11 <- matrix(NA, nrow = 19, ncol = 15)

for(Accuracy in c("high", "low")){

if(Accuracy == "high"){
  Se <- c(.95, .95)  # True assay sensitivities
  Sp <- c(.99, .99)  # True assay specificities
  i1 <- 0
  i2 <- 0
}

if(Accuracy == "low"){
  Se <- c(.85, .85)  # True assay sensitivities
  Sp <- c(.90, .90)  # True assay specificities
  i1 <- 13
  i2 <- 10
}

#========= Estimation with H2 protocol =========#

protocol <- H2
list.zmat <- lapply( rep(N, nsims), hier.alg.data, p=p,
                     design=protocol, Se=Se, Sp=Sp )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), p.pr=p.pr,
           Se1.pr=Se1.pr, Se2.pr=Se2.pr, Sp1.pr=Sp1.pr, Sp2.pr=Sp2.pr,
           postGit=G, method="Bayes", accuracy="unknown")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes, 
            p0=p0,delta0=delta0,N=N,S=length(protocol),p.pr=p.pr,
            Se1.pr=Se1.pr,Se2.pr=Se2.pr,Sp1.pr=Sp1.pr,Sp2.pr=Sp2.pr,
            emGit=G,emburn=burn,method="MAP",accuracy="unknown")

p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 4)
T <- rep(-9, nsims)

for(s in 1:nsims){
## POSTERIOR SAMPLING 
  res1 <- res.bayes[[s]]
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)
  delta.Mean[s, ] <- colMeans(res1$accuracy[-(1:burn),][pick,])
  delta.SE[s, ] <- apply(res1$accuracy[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence
  delta.MAP[s, ] <- res2$accuracy

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D10[i1+1:4, 1] <- round(colMeans(p.Mean),3)   # Est
Table.D10[i1+1:4, 2] <- round(apply(p.Mean,2,sd),4) # SD
Table.D10[i1+1:4, 3] <- round(colMeans(p.SE),4)     # SE

Table.D11[i2+1:4, 1] <- round(colMeans(delta.Mean),3)   # Est
Table.D11[i2+1:4, 2] <- round(apply(delta.Mean,2,sd),3) # SD
Table.D11[i2+1:4, 3] <- round(colMeans(delta.SE),3)     # SE

## MAP Estimation:
Table.D10[i1+6:9, 1] <- round(colMeans(p.MAP),3)    # Est
Table.D10[i1+6:9, 2] <- round(apply(p.MAP,2,sd),4)  # SD
Table.D10[i1+6:9, 3] <- round(colMeans(p.SE),4)     # SE

Table.D11[i2+6:9, 1] <- round(colMeans(delta.MAP),3)    # Est
Table.D11[i2+6:9, 2] <- round(apply(delta.MAP,2,sd),3)  # SD
Table.D11[i2+6:9, 3] <- round(colMeans(delta.SE),3)     # SE

## Avg. number of tests:
Table.D10[i1+11, 2] <- round( mean(T), 1 )    
Table.D10[i1+12, 2] <- round( sd(T), 1 )

#========= Estimation with H3 protocol =========#

protocol <- H3
list.zmat <- lapply( rep(N, nsims), hier.alg.data, p=p,
                     design=protocol, Se=Se, Sp=Sp )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), p.pr=p.pr,
           Se1.pr=Se1.pr, Se2.pr=Se2.pr, Sp1.pr=Sp1.pr, Sp2.pr=Sp2.pr,
           postGit=G, method="Bayes", accuracy="unknown")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes, 
            p0=p0,delta0=delta0,N=N,S=length(protocol),p.pr=p.pr,
            Se1.pr=Se1.pr,Se2.pr=Se2.pr,Sp1.pr=Sp1.pr,Sp2.pr=Sp2.pr,
            emGit=G,emburn=burn,method="MAP",accuracy="unknown")

p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 4)
T <- rep(-9, nsims)

for(s in 1:nsims){
## POSTERIOR SAMPLING 
  res1 <- res.bayes[[s]]
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)
  delta.Mean[s, ] <- colMeans(res1$accuracy[-(1:burn),][pick,])
  delta.SE[s, ] <- apply(res1$accuracy[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence
  delta.MAP[s, ] <- res2$accuracy

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D10[i1+1:4, 5] <- round(colMeans(p.Mean),3)   # Est
Table.D10[i1+1:4, 6] <- round(apply(p.Mean,2,sd),4) # SD
Table.D10[i1+1:4, 7] <- round(colMeans(p.SE),4)     # SE

Table.D11[i2+1:4, 5] <- round(colMeans(delta.Mean),3)   # Est
Table.D11[i2+1:4, 6] <- round(apply(delta.Mean,2,sd),3) # SD
Table.D11[i2+1:4, 7] <- round(colMeans(delta.SE),3)     # SE

## MAP Estimation:
Table.D10[i1+6:9, 5] <- round(colMeans(p.MAP),3)    # Est
Table.D10[i1+6:9, 6] <- round(apply(p.MAP,2,sd),4)  # SD
Table.D10[i1+6:9, 7] <- round(colMeans(p.SE),4)     # SE

Table.D11[i2+6:9, 5] <- round(colMeans(delta.MAP),3)    # Est
Table.D11[i2+6:9, 6] <- round(apply(delta.MAP,2,sd),3)  # SD
Table.D11[i2+6:9, 7] <- round(colMeans(delta.SE),3)     # SE

## Avg. number of tests:
Table.D10[i1+11, 6] <- round( mean(T), 1 )    
Table.D10[i1+12, 6] <- round( sd(T), 1 )

#========= Estimation with H4 protocol =========#

protocol <- H4
list.zmat <- lapply( rep(N, nsims), hier.alg.data, p=p,
                     design=protocol, Se=Se, Sp=Sp )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), p.pr=p.pr,
           Se1.pr=Se1.pr, Se2.pr=Se2.pr, Sp1.pr=Sp1.pr, Sp2.pr=Sp2.pr,
           postGit=G, method="Bayes", accuracy="unknown")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes, 
            p0=p0,delta0=delta0,N=N,S=length(protocol),p.pr=p.pr,
            Se1.pr=Se1.pr,Se2.pr=Se2.pr,Sp1.pr=Sp1.pr,Sp2.pr=Sp2.pr,
            emGit=G,emburn=burn,method="MAP",accuracy="unknown")

p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 4)
T <- rep(-9, nsims)

for(s in 1:nsims){
## POSTERIOR SAMPLING 
  res1 <- res.bayes[[s]]
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)
  delta.Mean[s, ] <- colMeans(res1$accuracy[-(1:burn),][pick,])
  delta.SE[s, ] <- apply(res1$accuracy[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence
  delta.MAP[s, ] <- res2$accuracy

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D10[i1+1:4, 9] <- round(colMeans(p.Mean),3)   # Est
Table.D10[i1+1:4, 10] <- round(apply(p.Mean,2,sd),4) # SD
Table.D10[i1+1:4, 11] <- round(colMeans(p.SE),4)     # SE

Table.D11[i2+1:4, 9] <- round(colMeans(delta.Mean),3)   # Est
Table.D11[i2+1:4, 10] <- round(apply(delta.Mean,2,sd),3) # SD
Table.D11[i2+1:4, 11] <- round(colMeans(delta.SE),3)     # SE

## MAP Estimation:
Table.D10[i1+6:9, 9] <- round(colMeans(p.MAP),3)    # Est
Table.D10[i1+6:9, 10] <- round(apply(p.MAP,2,sd),4)  # SD
Table.D10[i1+6:9, 11] <- round(colMeans(p.SE),4)     # SE

Table.D11[i2+6:9, 9] <- round(colMeans(delta.MAP),3)    # Est
Table.D11[i2+6:9, 10] <- round(apply(delta.MAP,2,sd),3)  # SD
Table.D11[i2+6:9, 11] <- round(colMeans(delta.SE),3)     # SE

## Avg. number of tests:
Table.D10[i1+11, 10] <- round( mean(T), 1 )    
Table.D10[i1+12, 10] <- round( sd(T), 1 )

#========= Estimation with AT protocol =========#

protocol <- AT
list.zmat <- lapply( rep(N, nsims), array.2dim.data, p=p,
                     design=protocol, Se=Se, Sp=Sp )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), p.pr=p.pr,
           Se1.pr=Se1.pr, Se2.pr=Se2.pr, Sp1.pr=Sp1.pr, Sp2.pr=Sp2.pr,
           postGit=G, method="Bayes", accuracy="unknown")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes, 
            p0=p0,delta0=delta0,N=N,S=length(protocol),p.pr=p.pr,
            Se1.pr=Se1.pr,Se2.pr=Se2.pr,Sp1.pr=Sp1.pr,Sp2.pr=Sp2.pr,
            emGit=G,emburn=burn,method="MAP",accuracy="unknown")

p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 4)
T <- rep(-9, nsims)

for(s in 1:nsims){
## POSTERIOR SAMPLING 
  res1 <- res.bayes[[s]]
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)
  delta.Mean[s, ] <- colMeans(res1$accuracy[-(1:burn),][pick,])
  delta.SE[s, ] <- apply(res1$accuracy[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence
  delta.MAP[s, ] <- res2$accuracy

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.D10[i1+1:4, 13] <- round(colMeans(p.Mean),3)   # Est
Table.D10[i1+1:4, 14] <- round(apply(p.Mean,2,sd),4) # SD
Table.D10[i1+1:4, 15] <- round(colMeans(p.SE),4)     # SE

Table.D11[i2+1:4, 13] <- round(colMeans(delta.Mean),3)   # Est
Table.D11[i2+1:4, 14] <- round(apply(delta.Mean,2,sd),3) # SD
Table.D11[i2+1:4, 15] <- round(colMeans(delta.SE),3)     # SE

## MAP Estimation:
Table.D10[i1+6:9, 13] <- round(colMeans(p.MAP),3)    # Est
Table.D10[i1+6:9, 14] <- round(apply(p.MAP,2,sd),4)  # SD
Table.D10[i1+6:9, 15] <- round(colMeans(p.SE),4)     # SE

Table.D11[i2+6:9, 13] <- round(colMeans(delta.MAP),3)    # Est
Table.D11[i2+6:9, 14] <- round(apply(delta.MAP,2,sd),3)  # SD
Table.D11[i2+6:9, 15] <- round(colMeans(delta.SE),3)     # SE

## Avg. number of tests:
Table.D10[i1+11, 14] <- round( mean(T), 1 )    
Table.D10[i1+12, 14] <- round( sd(T), 1 )

}

stopCluster( cl )

## Preparing & saving Table D10
dfr1 <- data.frame( cbind(
     c("","","","","","High accuracy","","","","","","", 
	   "",  "","","","","","Low accuracy","","","","","",""),
     c("", "", "", "p00=0.95", "p10=0.02", "p01=0.02", "p11=0.01", "", "", "", "", "", "", 
	   "", "", "", "p00=0.95", "p10=0.02", "p01=0.02", "p11=0.01", "", "", "", "", ""),  
     c("Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T", "", 
	   "Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T")))

dfr1 <- data.frame( dfr1, Table.D10 )
dfr1 <- sapply(dfr1, as.character)
dfr1[is.na(dfr1)] <- " "
Table_D10 <- as.data.frame(dfr1)
colnames( Table_D10 ) <- c("", "","", "Est", "SD", "SE", "", "Est", "SD", "SE", 
                           "", "Est", "SD", "SE", "", "Est", "SD", "SE")				   
saveRDS(Table_D10, "./simulation/intermediate_results/Table_D10.rds")

## Preparing & saving Table D11
dfr2 <- data.frame( cbind(
     c("","","","","","High accuracy","","","","","","", 
	   "","","","Low accuracy","","",""),
     c("", "", "", "Se:(1)1=0.95", "Se:(1)2=0.95", "Sp:(1)1=0.99", "Sp:(1)2=0.99", "", "", "",  
	   "", "", "", "Se:(1)1=0.85", "Se:(1)2=0.85", "Sp:(1)1=0.90", "Sp:(1)2=0.90", "", ""),  
     c("Mean", "", "", "", "", "MAP", "", "", "", "",  
	   "Mean", "", "", "", "", "MAP", "", "", "")))

dfr2 <- data.frame( dfr2, Table.D11 )
dfr2 <- sapply(dfr2, as.character)
dfr2[is.na(dfr2)] <- " "

Table_D11 <- as.data.frame(dfr2)
colnames( Table_D11 ) <- c("", "","", "Est", "SD", "SE", "", "Est", "SD", "SE", 
                           "", "Est", "SD", "SE", "", "Est", "SD", "SE")
						   
saveRDS(Table_D11, "./simulation/intermediate_results/Table_D11.rds")

