
##################################################

# 'ncores' is the number of cores to be used
# for parallel computing in Windows with 64-bit R:
ncores <- 10

# 'nsims' is the number of simulated data sets:
nsims <- 500

## Importing source code & loading packages:
library(MCMCpack)
library("parallel")
library(Rcpp)
library(RcppArmadillo)
library(binGroup2)

source("./source_code/MainFunctions.R")
source("./source_code/SupportPrograms.R")
source("./source_code/GT_sim_fns.R")
sourceCpp("./source_code/Hierarchical.cpp")

# Making clusters for parallel computing:
cl <- makeCluster(ncores)
clusterExport(cl, "hier.alg.data")
clusterExport(cl, "array.2dim.data")
clusterExport(cl, "PostUnknownAssayAcr")
clusterExport(cl, "poolMembTracker")
clusterExport(cl, "EmUnknownAssayAcr")
clusterExport(cl, "EmultGibbs2")

##################################################

## This block calculates the optimal pool sizes
## shown in Table 1.

S <- 4
se.mat <- matrix(0.95, S, 2)
sp.mat <- matrix(0.99, S, 2)
opt.psz <- array(NA, c(S, S, 2), )

for( d in 1:2 ){
  if(d == 1){ 
    p <- c(.95, .02, .02, .01)     ## Config-I
  }
  if( d == 2){
    p <- c(.990, .004, .004, .002) ## Config-II
  }

  ## Hierarchical testing
  for(s in 2:S){
    designs <- read.table( paste("./source_code/", "psz_s", s, ".csv", sep="") )
    Se <- se.mat[1:s, ]
    Sp <- sp.mat[1:s, ]
    hier.et <- rep(-9, nrow(designs))
    for(i in 1:nrow(designs)){
      hier.et[i] <- EFF(p=p, S=s, SE=Se, SP=Sp, 
                        ns=as.numeric(designs[i, ]))
    }
    opt.psz[s-1, 1:s, d] <- as.numeric( designs[hier.et == min(hier.et), ] )
  }

  ## Array testing
  max.ary.sz <- 30     ## max array size
  designs <- 2:max.ary.sz
  at.et <- rep(-9, length(designs))
  for(k in designs){
    res <- OTC2(algorithm = "A2", p.vec = p, Se = Se[1:2, ], 
                Sp = Sp[1:2, ], group.sz = k, trace = FALSE, print.time = FALSE)
    at.et[k-1] <- res[5]$opt.ET$value
  }
  opt.psz[S, 1:2, d] <- designs[at.et == min(at.et)]
}

Table_1 <- data.frame( c("H2", "H3", "H4", "AT"), opt.psz[ , ,1], 
                       c("", "", "", ""), 
                       c("H2", "H3", "H4", "AT"), opt.psz[ , ,2])
colnames( Table_1 ) <- c("Protocol", "Pool sizes", "", "", "", "", "Protocol", "Pool sizes", "", "", "")
Table_1[ is.na(Table_1) ] <- ""

saveRDS(Table_1, "./simulation/intermediate_results/Table_1.rds")

##################################################

## The code below reproduces the simulation 
## results depicted in Tables 2 and 3.

# Setting seed:
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(123)
clusterSetRNGStream(cl = cl, iseed = 123)

# Simulation parameters:
N <- 5000         # Sample size
Se <- c(.95, .95)  # True sensitivities for diseases 1 & 2
Sp <- c(.99, .99)  # True specificities for diseases 1 & 2

## Specify the prior hyperparameters as shown below. 
## Flat priors have been used for all parameters.
p.pr <- rep(1, 4)
Se1.pr <- c(1, 1)
Se2.pr <- c(1, 1)
Sp1.pr <- c(1, 1)
Sp2.pr <- c(1, 1)

G <- 12000               # number of Gibbs iterates
burn <- 2000             # burn-in period
pick <- seq(1, 10000, 5)   # keeping every 5th

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
Table.2 <- matrix(NA, nrow = 25, ncol = 15)
Table.3 <- matrix(NA, nrow = 19, ncol = 15)

for(Configuration in c("I", "II")){

if(Configuration == "I"){
  ## Configuration I
  p <- c(.95, .02, .02, .01)     # True value of p
  H2 <- c(5, 1)
  H3 <- c(9, 3, 1)
  H4 <- c(18, 6, 3, 1)
  AT <- c(11, 11, 1)  # For AT, c(11,11,1), not c(11,11)
                    # because tests are done in 3 stages:
                    # row pool, col pool, & individual testing
  i1 <- 0
  i2 <- 0
}

if(Configuration == "II"){
  ## Configuration II
  p <- c(.990, .004, .004, .002) # True value of p
  H2 <- c(11, 1)
  H3 <- c(25, 5, 1)
  H4 <- c(48, 12, 4, 1)
  AT <- c(29, 29, 1)
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
            p0=p0, delta0=delta0, N=N, S=length(protocol), p.pr=p.pr,
            Se1.pr=Se1.pr, Se2.pr=Se2.pr, Sp1.pr=Sp1.pr, Sp2.pr=Sp2.pr,
            emGit=G, emburn=burn, method="MAP", accuracy="unknown")

p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 4)
T <- rep(-9, nsims)

for(s in 1:nsims){
## POSTERIOR SAMPLING 
  res1 <- res.bayes[[s]]
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn), ][pick,])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn), ][pick, ], 2, sd)
  delta.Mean[s, ] <- colMeans(res1$accuracy[-(1:burn), ][pick, ])
  delta.SE[s, ] <- apply(res1$accuracy[-(1:burn), ][pick, ], 2, sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence
  delta.MAP[s, ] <- res2$accuracy

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.2[i1+1:4, 1] <- round(colMeans(p.Mean), 3)   # Est
Table.2[i1+1:4, 2] <- round(apply(p.Mean, 2, sd), 4) # SD
Table.2[i1+1:4, 3] <- round(colMeans(p.SE), 4)     # SE

Table.3[i2+1:4, 1] <- round(colMeans(delta.Mean), 3)   # Est
Table.3[i2+1:4, 2] <- round(apply(delta.Mean, 2, sd), 3) # SD
Table.3[i2+1:4, 3] <- round(colMeans(delta.SE), 3)     # SE

## MAP Estimation:
Table.2[i1+6:9, 1] <- round(colMeans(p.MAP), 3)    # Est
Table.2[i1+6:9, 2] <- round(apply(p.MAP, 2, sd), 4)  # SD
Table.2[i1+6:9, 3] <- round(colMeans(p.SE), 4)     # SE

Table.3[i2+6:9, 1] <- round(colMeans(delta.MAP), 3)    # Est
Table.3[i2+6:9, 2] <- round(apply(delta.MAP, 2, sd), 3)  # SD
Table.3[i2+6:9, 3] <- round(colMeans(delta.SE), 3)     # SE

## Avg. number of tests:
Table.2[i1+11, 1] <- round( mean(T), 1 )    
Table.2[i1+12, 1] <- round( sd(T), 1 )

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
            p0=p0, delta0=delta0, N=N, S=length(protocol), p.pr=p.pr,
            Se1.pr=Se1.pr, Se2.pr=Se2.pr, Sp1.pr=Sp1.pr, Sp2.pr=Sp2.pr,
            emGit=G, emburn=burn, method="MAP", accuracy="unknown")

p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 4)
T <- rep(-9, nsims)

for(s in 1:nsims){
## POSTERIOR SAMPLING 
  res1 <- res.bayes[[s]]
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn), ][pick,])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn), ][pick, ], 2, sd)
  delta.Mean[s, ] <- colMeans(res1$accuracy[-(1:burn), ][pick, ])
  delta.SE[s, ] <- apply(res1$accuracy[-(1:burn), ][pick, ], 2, sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence
  delta.MAP[s, ] <- res2$accuracy

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.2[i1+1:4, 5] <- round(colMeans(p.Mean), 3)   # Est
Table.2[i1+1:4, 6] <- round(apply(p.Mean, 2, sd), 4) # SD
Table.2[i1+1:4, 7] <- round(colMeans(p.SE), 4)     # SE

Table.3[i2+1:4, 5] <- round(colMeans(delta.Mean), 3)   # Est
Table.3[i2+1:4, 6] <- round(apply(delta.Mean,2, sd), 3) # SD
Table.3[i2+1:4, 7] <- round(colMeans(delta.SE), 3)     # SE

## MAP Estimation:
Table.2[i1+6:9, 5] <- round(colMeans(p.MAP), 3)    # Est
Table.2[i1+6:9, 6] <- round(apply(p.MAP, 2, sd), 4)  # SD
Table.2[i1+6:9, 7] <- round(colMeans(p.SE), 4)     # SE

Table.3[i2+6:9, 5] <- round(colMeans(delta.MAP), 3)    # Est
Table.3[i2+6:9, 6] <- round(apply(delta.MAP, 2, sd), 3)  # SD
Table.3[i2+6:9, 7] <- round(colMeans(delta.SE), 3)     # SE

## Avg. number of tests:
Table.2[i1+11, 5] <- round( mean(T), 1 )    
Table.2[i1+12, 5] <- round( sd(T), 1 )

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
            p0=p0, delta0=delta0, N=N, S=length(protocol), p.pr=p.pr,
            Se1.pr=Se1.pr, Se2.pr=Se2.pr, Sp1.pr=Sp1.pr, Sp2.pr=Sp2.pr,
            emGit=G, emburn=burn, method="MAP", accuracy="unknown")

p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 4)
T <- rep(-9, nsims)

for(s in 1:nsims){
## POSTERIOR SAMPLING 
  res1 <- res.bayes[[s]]
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn), ][pick, ])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn), ][pick, ], 2, sd)
  delta.Mean[s, ] <- colMeans(res1$accuracy[-(1:burn), ][pick, ])
  delta.SE[s, ] <- apply(res1$accuracy[-(1:burn), ][pick, ], 2, sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence
  delta.MAP[s, ] <- res2$accuracy

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.2[i1+1:4, 9] <- round(colMeans(p.Mean), 3)   # Est
Table.2[i1+1:4, 10] <- round(apply(p.Mean, 2, sd), 4) # SD
Table.2[i1+1:4, 11] <- round(colMeans(p.SE), 4)     # SE

Table.3[i2+1:4, 9] <- round(colMeans(delta.Mean), 3)   # Est
Table.3[i2+1:4, 10] <- round(apply(delta.Mean, 2, sd), 3) # SD
Table.3[i2+1:4, 11] <- round(colMeans(delta.SE), 3)     # SE

## MAP Estimation:
Table.2[i1+6:9, 9] <- round(colMeans(p.MAP), 3)    # Est
Table.2[i1+6:9, 10] <- round(apply(p.MAP, 2, sd), 4)  # SD
Table.2[i1+6:9, 11] <- round(colMeans(p.SE), 4)     # SE

Table.3[i2+6:9, 9] <- round(colMeans(delta.MAP), 3)    # Est
Table.3[i2+6:9, 10] <- round(apply(delta.MAP, 2, sd), 3)  # SD
Table.3[i2+6:9, 11] <- round(colMeans(delta.SE), 3)     # SE

## Avg. number of tests:
Table.2[i1+11, 9] <- round( mean(T), 1 )    
Table.2[i1+12, 9] <- round( sd(T), 1 )

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
            p0=p0, delta0=delta0, N=N, S=length(protocol), p.pr=p.pr,
            Se1.pr=Se1.pr, Se2.pr=Se2.pr, Sp1.pr=Sp1.pr, Sp2.pr=Sp2.pr,
            emGit=G, emburn=burn, method="MAP", accuracy="unknown")

p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 4)
T <- rep(-9, nsims)

for(s in 1:nsims){
## POSTERIOR SAMPLING 
  res1 <- res.bayes[[s]]
  p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn), ][pick, ])
  p.SE[s, ] <- apply(res1$prevalence[-(1:burn), ][pick, ], 2, sd)
  delta.Mean[s, ] <- colMeans(res1$accuracy[-(1:burn), ][pick, ])
  delta.SE[s, ] <- apply(res1$accuracy[-(1:burn), ][pick, ], 2, sd)

## MAP ESTIMATION
  res2 <- res.map[[s]]
  p.MAP[s, ] <- res2$prevalence
  delta.MAP[s, ] <- res2$accuracy

## Number of tests
  T[s] <- nrow( list.zmat[[s]] )
}

## Mean Estimation:
Table.2[i1+1:4, 13] <- round(colMeans(p.Mean), 3)   # Est
Table.2[i1+1:4, 14] <- round(apply(p.Mean, 2, sd), 4) # SD
Table.2[i1+1:4, 15] <- round(colMeans(p.SE), 4)     # SE

Table.3[i2+1:4, 13] <- round(colMeans(delta.Mean), 3)   # Est
Table.3[i2+1:4, 14] <- round(apply(delta.Mean, 2, sd), 3) # SD
Table.3[i2+1:4, 15] <- round(colMeans(delta.SE), 3)     # SE

## MAP Estimation:
Table.2[i1+6:9, 13] <- round(colMeans(p.MAP), 3)    # Est
Table.2[i1+6:9, 14] <- round(apply(p.MAP, 2, sd), 4)  # SD
Table.2[i1+6:9, 15] <- round(colMeans(p.SE), 4)     # SE

Table.3[i2+6:9, 13] <- round(colMeans(delta.MAP), 3)    # Est
Table.3[i2+6:9, 14] <- round(apply(delta.MAP, 2, sd), 3)  # SD
Table.3[i2+6:9, 15] <- round(colMeans(delta.SE), 3)     # SE

## Avg. number of tests:
Table.2[i1+11, 13] <- round( mean(T), 1 )    
Table.2[i1+12, 13] <- round( sd(T), 1 )

}

stopCluster( cl )

## Preparing & saving Table 2
dfr1 <- data.frame( cbind(
     c("","","","","","Config-I","","","","","","", 
	   "",  "","","","","","Config-II","","","","","",""),
     c("", "", "", "p00=0.95", "p10=0.02", "p01=0.02", "p11=0.01", "", "", "", "", "", "", 
	   "", "", "", "p00=0.990", "p10=0.004", "p01=0.004", "p11=0.002", "", "", "", "", ""),  
     c("Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T", "", 
	   "Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T")))

dfr1 <- data.frame( dfr1, Table.2 )
dfr1 <- sapply(dfr1, as.character)
dfr1[is.na(dfr1)] <- " "
Table_2 <- as.data.frame(dfr1)
colnames( Table_2 ) <- c("", "","", "Est", "SD", "SE", "", "Est", "SD", "SE", 
                           "", "Est", "SD", "SE", "", "Est", "SD", "SE")				   

## Mean
T <- as.numeric( c( Table_2[11, 4], Table_2[11, 8], Table_2[11, 12], Table_2[11, 16] ) )
redc <- round( 100*(T[1] - T[-1])/T[1], 1 )
T.redc <- c( T[1],
             paste( T[2], "(", redc[1], "%", ")", sep="" ),
             paste( T[3], "(", redc[2], "%", ")", sep="" ),
             paste( T[4], "(", redc[3], "%", ")", sep="" )
            )
Table_2[11, 8] <- T.redc[2]
Table_2[11, 12] <- T.redc[3]
Table_2[11, 16] <- T.redc[4]

## MAP
T <- as.numeric( c( Table_2[24, 4], Table_2[24, 8], Table_2[24, 12], Table_2[24, 16] ) )
redc <- round( 100*(T[1] - T[-1])/T[1], 1 )
T.redc <- c( T[1],
             paste( T[2], "(", redc[1], "%", ")", sep="" ),
             paste( T[3], "(", redc[2], "%", ")", sep="" ),
             paste( T[4], "(", redc[3], "%", ")", sep="" )
            )
Table_2[24, 8] <- T.redc[2]
Table_2[24, 12] <- T.redc[3]
Table_2[24, 16] <- T.redc[4]
saveRDS(Table_2, "./simulation/intermediate_results/Table_2.rds")


## Preparing & saving Table 3
dfr2 <- data.frame( cbind(
     c("","","","","","Config-I","","","","","","", 
	   "","","","Config-II","","",""),
     c("", "", "", "Se:(1)1=0.95", "Se:(1)2=0.95", "Sp:(1)1=0.99", "Sp:(1)2=0.99", "", "", "",  
	   "", "", "", "Se:(1)1=0.95", "Se:(1)2=0.95", "Sp:(1)1=0.99", "Sp:(1)2=0.99", "", ""),  
     c("Mean", "", "", "", "", "MAP", "", "", "", "",  
	   "Mean", "", "", "", "", "MAP", "", "", "")))

dfr2 <- data.frame( dfr2, Table.3 )
dfr2 <- sapply(dfr2, as.character)
dfr2[is.na(dfr2)] <- " "

Table_3 <- as.data.frame(dfr2)
colnames( Table_3 ) <- c("", "","", "Est", "SD", "SE", "", "Est", "SD", "SE", 
                           "", "Est", "SD", "SE", "", "Est", "SD", "SE")
						   
saveRDS(Table_3, "./simulation/intermediate_results/Table_3.rds")

