
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
clusterExport(cl, "PostUnknownAssayAcr_L2")
clusterExport(cl, "poolMembTracker2_L2")
clusterExport(cl, "EmUnknownAssayAcr_L2")
clusterExport(cl, "EmultGibbs2_L2")


##################################################

## This block calculates the optimal pool sizes
## shown in Table D.2 of the Supporting Information.

S <- 4

# For H2
Se.H2 <- rbind(c(.95,.95), c(.98,.98) )
Sp.H2 <- rbind(c(.98,.98), c(.99,.99) )
			  
# For H3
Se.H3 <- rbind(c(.95,.95), c(.95,.95), c(.98,.98) )
Sp.H3 <- rbind(c(.98,.98), c(.98,.98), c(.99,.99) )

# For H4
Se.H4 <- rbind(c(.95,.95), c(.95,.95), c(.95,.95), c(.98,.98) )
Sp.H4 <- rbind(c(.98,.98), c(.98,.98), c(.98,.98), c(.99,.99) )

se.list <- list( Se.H2, Se.H3, Se.H4 )
sp.list <- list( Sp.H2, Sp.H3, Sp.H4 )

# For AT
Se.AT <- rbind(c(.95,.95), c(.98,.98) )
Sp.AT <- rbind(c(.98,.98), c(.99,.99) )

opt.psz <- array(NA, c(S, S, 2), )

for( d in 1:2 ){
  if(d == 1){ 
    p <- c(.95,.02,.02,.01)     ## Config-I
  }
  if( d == 2){
    p <- c(.990,.004,.004,.002) ## Config-II
  }

  ## Hierarchical testing
  for(s in 2:S){
    designs <- read.table( paste("./source_code/","psz_s", s, ".csv", sep="") )
    Se <- se.list[[ s-1 ]]
    Sp <- sp.list[[ s-1 ]]
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
    res <- OTC2(algorithm = "A2", p.vec = p, Se = t( Se.AT ), 
                Sp = t( Sp.AT ), group.sz = k, trace = FALSE, print.time = FALSE)
    at.et[k-1] <- res[5]$opt.ET$value
  }
  opt.psz[S, 1:2, d] <- designs[at.et == min(at.et)]
}

Table_D2 <- data.frame( c("H2", "H3", "H4", "AT"), opt.psz[ , ,1], 
                       c("","","",""), 
                       c("H2", "H3", "H4", "AT"), opt.psz[ , ,2])
colnames( Table_D2 ) <- c("Protocol", "Pool sizes", "", "", "", 
                          "", "Protocol", "Pool sizes", "", "", "")
Table_D2[ is.na(Table_D2) ] <- ""

saveRDS(Table_D2, "./simulation/intermediate_results/Table_D2.rds")

##################################################

# Sample size
N <- 5000         # Sample size

# Setting seed:
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed(123)
clusterSetRNGStream(cl = cl, iseed = 123)

G <- 12000               # number of Gibbs iterates
burn <- 2000             # burn-in period
pick <- seq(1,10000,5)   # keeping every 5th

## Convergence tolerance for the EM algorithm
epsilon <- 0.001

## Maximum number of iterations for the EM algorithm. 
emmaxit <- 200

options( scipen = 999 )

## Remarks about the initial values:
## --------------------------------
## The proposed estimation techniques require an initial value of 
## p and delta. These values can be specified from historical, pilot
## study, or any contextual estimate. For p, we use p0=c(.92,.05,.02,.01). 
## For all sensitivities and specificities, we use .96 and .98,
## respectively, as done with other simulations in this paper.
p0 <- c(.92, .05, .02, .01)
delta0 <- c(.96, .96, .98, .98, .96, .96, .98, .98)

## Note: We perform two simulations. The first one 
## reproduces the results shown in Tables D.3-D.5, while 
## the second reproduces Tables D.6-D.8 results.


#### SIMULATION I
#### Flat priors for all parameters (Tables D.3-D.5):

p.pr <- rep(1, 4)
se.pr <- matrix(1, 2, 4)
sp.pr <- matrix(1, 2, 4)

Table.D3 <- matrix(NA, nrow = 25, ncol = 15)
Table.D4 <- matrix(NA, nrow = 17, ncol = 15)
Table.D5 <- matrix(NA, nrow = 17, ncol = 15)

for(Configuration in c("I", "II")){

if(Configuration == "I"){
  ## Configuration I
  p <- c(.95,.02,.02,.01)     # True value of p
  H2 <- c(5, 1)
  H3 <- c(9, 3, 1)
  H4 <- c(18, 6, 3, 1)
  AT <- c(11, 11, 1)
  i1 <- 0
  i2 <- 0
}

if(Configuration == "II"){
  ## Configuration II
  p <- c(.990,.004,.004,.002) # True value of p
  H2 <- c(11, 1)
  H3 <- c(25, 5, 1)
  H4 <- c(36, 12, 4, 1)
  AT <- c(28, 28, 1)
  i1 <- 13
  i2 <- 10
}

#========= Estimation with H2 protocol =========#

protocol <- H2
list.zmat <- lapply( rep(N, nsims), hier.alg.data_L2, p=p,
                     design=protocol, Se=Se.H2, Sp=Sp.H2 )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           postGit=G, method="Bayes")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           emGit=G, emburn=burn, method="MAP")


p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 8)
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

## Mean:
Table.D3[i1+1:4, 1] <- round(colMeans(p.Mean), 3)   # Est
Table.D3[i1+1:4, 2] <- round(apply(p.Mean,2,sd), 4) # SD
Table.D3[i1+1:4, 3] <- round(colMeans(p.SE), 4)     # SE

## MAP:
Table.D3[i1+6:9, 1] <- round(colMeans(p.MAP), 3)    # Est
Table.D3[i1+6:9, 2] <- round(apply(p.MAP,2,sd), 4)  # SD
Table.D3[i1+6:9, 3] <- round(colMeans(p.SE), 4)     # SE

if( Configuration == "I" ){
  # Mean
  Table.D4[1:8, 1] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D4[1:8, 2] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D4[1:8, 3] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D4[10:17, 1] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D4[10:17, 2] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D4[10:17, 3] <- round(colMeans(delta.SE), 3)     # SE
}

if(Configuration == "II"){
  # Mean
  Table.D5[1:8, 1] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D5[1:8, 2] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D5[1:8, 3] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D5[10:17, 1] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D5[10:17, 2] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D5[10:17, 3] <- round(colMeans(delta.SE), 3)     # SE
}

## Avg. number of tests:
Table.D3[i1+11, 1] <- round( mean(T), 1 )    
Table.D3[i1+12, 1] <- round( sd(T), 1 )


#========= Estimation with H3 protocol =========#

protocol <- H3
list.zmat <- lapply( rep(N, nsims), hier.alg.data_L2, p=p,
                     design=protocol, Se=Se.H3, Sp=Sp.H3 )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           postGit=G, method="Bayes")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           emGit=G, emburn=burn, method="MAP")


p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 8)
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

## Mean:
Table.D3[i1+1:4, 5] <- round(colMeans(p.Mean), 3)   # Est
Table.D3[i1+1:4, 6] <- round(apply(p.Mean,2,sd), 4) # SD
Table.D3[i1+1:4, 7] <- round(colMeans(p.SE), 4)     # SE

## MAP:
Table.D3[i1+6:9, 5] <- round(colMeans(p.MAP), 3)    # Est
Table.D3[i1+6:9, 6] <- round(apply(p.MAP,2,sd), 4)  # SD
Table.D3[i1+6:9, 7] <- round(colMeans(p.SE), 4)     # SE

if( Configuration == "I" ){
  # Mean
  Table.D4[1:8, 5] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D4[1:8, 6] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D4[1:8, 7] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D4[10:17, 5] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D4[10:17, 6] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D4[10:17, 7] <- round(colMeans(delta.SE), 3)     # SE
}

if(Configuration == "II"){
  # Mean
  Table.D5[1:8, 5] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D5[1:8, 6] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D5[1:8, 7] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D5[10:17, 5] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D5[10:17, 6] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D5[10:17, 7] <- round(colMeans(delta.SE), 3)     # SE
}

## Avg. number of tests:
Table.D3[i1+11, 5] <- round( mean(T), 1 )    
Table.D3[i1+12, 5] <- round( sd(T), 1 )


#========= Estimation with H4 protocol =========#

protocol <- H4
list.zmat <- lapply( rep(N, nsims), hier.alg.data_L2, p=p,
                     design=protocol, Se=Se.H4, Sp=Sp.H4 )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           postGit=G, method="Bayes")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           emGit=G, emburn=burn, method="MAP")


p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 8)
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

## Mean:
Table.D3[i1+1:4, 9] <- round(colMeans(p.Mean), 3)   # Est
Table.D3[i1+1:4, 10] <- round(apply(p.Mean,2,sd), 4) # SD
Table.D3[i1+1:4, 11] <- round(colMeans(p.SE), 4)     # SE

## MAP:
Table.D3[i1+6:9, 9] <- round(colMeans(p.MAP), 3)    # Est
Table.D3[i1+6:9, 10] <- round(apply(p.MAP,2,sd), 4)  # SD
Table.D3[i1+6:9, 11] <- round(colMeans(p.SE), 4)     # SE

if( Configuration == "I" ){
  # Mean
  Table.D4[1:8, 9] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D4[1:8, 10] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D4[1:8, 11] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D4[10:17, 9] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D4[10:17, 10] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D4[10:17, 11] <- round(colMeans(delta.SE), 3)     # SE
}

if(Configuration == "II"){
  # Mean
  Table.D5[1:8, 9] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D5[1:8, 10] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D5[1:8, 11] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D5[10:17, 9] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D5[10:17, 10] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D5[10:17, 11] <- round(colMeans(delta.SE), 3)     # SE
}

## Avg. number of tests:
Table.D3[i1+11, 9] <- round( mean(T), 1 )    
Table.D3[i1+12, 9] <- round( sd(T), 1 )


#========= Estimation with AT protocol =========#

protocol <- AT
list.zmat <- lapply( rep(N, nsims), array.2dim.data_L2, p=p,
                     design=protocol, Se=Se.AT, Sp=Sp.AT )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           postGit=G, method="Bayes")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           emGit=G, emburn=burn, method="MAP")


p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 8)
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

## Mean:
Table.D3[i1+1:4, 13] <- round(colMeans(p.Mean), 3)   # Est
Table.D3[i1+1:4, 14] <- round(apply(p.Mean,2,sd), 4) # SD
Table.D3[i1+1:4, 15] <- round(colMeans(p.SE), 4)     # SE

## MAP:
Table.D3[i1+6:9, 13] <- round(colMeans(p.MAP), 3)    # Est
Table.D3[i1+6:9, 14] <- round(apply(p.MAP,2,sd), 4)  # SD
Table.D3[i1+6:9, 15] <- round(colMeans(p.SE), 4)     # SE

if( Configuration == "I" ){
  # Mean
  Table.D4[1:8, 13] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D4[1:8, 14] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D4[1:8, 15] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D4[10:17, 13] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D4[10:17, 14] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D4[10:17, 15] <- round(colMeans(delta.SE), 3)     # SE
}

if(Configuration == "II"){
  # Mean
  Table.D5[1:8, 13] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D5[1:8, 14] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D5[1:8, 15] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D5[10:17, 13] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D5[10:17, 14] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D5[10:17, 15] <- round(colMeans(delta.SE), 3)     # SE
}

## Avg. number of tests:
Table.D3[i1+11, 13] <- round( mean(T), 1 )    
Table.D3[i1+12, 13] <- round( sd(T), 1 )

}

## Preparing & saving Table D.3
dfr1 <- data.frame( cbind(
     c("","","","","","Config-I","","","","","","", 
	   "",  "","","","","","Config-II","","","","","",""),
     c("", "", "", "p00=0.95", "p10=0.02", "p01=0.02", "p11=0.01", "", "", "", "", "", "", 
	   "", "", "", "p00=0.990", "p10=0.004", "p01=0.004", "p11=0.002", "", "", "", "", ""),  
     c("Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T", "", 
	   "Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T")))

dfr1 <- data.frame( dfr1, Table.D3 )
dfr1 <- sapply(dfr1, as.character)
dfr1[is.na(dfr1)] <- " "
Table_D3 <- as.data.frame(dfr1)
colnames( Table_D3 ) <- c("", "", "", "Est", "SD", "SE", "", "Est", "SD", "SE", 
                           "", "Est", "SD", "SE", "", "Est", "SD", "SE")				   

saveRDS(Table_D3, "./simulation/intermediate_results/Table_D3.rds")


## Preparing & saving Table D.4 & D.5
dfr2 <- data.frame( cbind(
     c("", "", "", "", "Se:(1)1=0.95", "Se:(1)2=0.95", "Sp:(1)1=0.98", "Sp:(1)2=0.98",
       "Se:(2)1=0.98", "Se:(2)2=0.98", "Sp:(2)1=0.99", "Sp:(2)2=0.99", "", "", "", "", ""),
     c("", "", "", "Mean", "", "", "", "", "", "", "", "", "", "MAP", "", "", "")))

dfr3 <- data.frame( dfr2, Table.D4 )
dfr3 <- sapply(dfr3, as.character)
dfr3[is.na(dfr3)] <- " "

Table_D4 <- as.data.frame(dfr3)
colnames( Table_D4 ) <- c("", "", "Est", "SD", "SE", "", "Est", "SD", "SE", 
                           "", "Est", "SD", "SE", "", "Est", "SD", "SE")
						   
saveRDS(Table_D4, "./simulation/intermediate_results/Table_D4.rds")


dfr4 <- data.frame( dfr2, Table.D5 )
dfr4 <- sapply(dfr4, as.character)
dfr4[is.na(dfr4)] <- " "

Table_D5 <- as.data.frame(dfr4)
colnames( Table_D5 ) <- c("", "", "Est", "SD", "SE", "", "Est", "SD", "SE", 
                           "", "Est", "SD", "SE", "", "Est", "SD", "SE")
						   
saveRDS(Table_D5, "./simulation/intermediate_results/Table_D5.rds")



#### SIMULATION II
#### Flat prior for p but informative for Se & Sp (Tables D.6-D.8):

p.pr <- rep(1, 4)
se.pr <- rbind(c(109.0, 6.7, 109.0, 6.7),  # for stage 1
               c(100.0, 3.0, 100.0, 3.0))  # for stage 2

sp.pr <- rbind(c(100.0, 3.0, 100.0, 3.0),  # for stage 1
               c( 55.2, 1.6,  55.2, 1.6))  # for stage 2

Table.D6 <- matrix(NA, nrow = 25, ncol = 15)
Table.D7 <- matrix(NA, nrow = 17, ncol = 15)
Table.D8 <- matrix(NA, nrow = 17, ncol = 15)

for(Configuration in c("I", "II")){

if(Configuration == "I"){
  ## Configuration I
  p <- c(.95,.02,.02,.01)     # True value of p
  H2 <- c(5, 1)
  H3 <- c(9, 3, 1)
  H4 <- c(18, 6, 3, 1)
  AT <- c(11, 11, 1)
  i1 <- 0
  i2 <- 0
}

if(Configuration == "II"){
  ## Configuration II
  p <- c(.990,.004,.004,.002) # True value of p
  H2 <- c(11, 1)
  H3 <- c(25, 5, 1)
  H4 <- c(36, 12, 4, 1)
  AT <- c(28, 28, 1)
  i1 <- 13
  i2 <- 10
}

#========= Estimation with H2 protocol =========#

protocol <- H2
list.zmat <- lapply( rep(N, nsims), hier.alg.data_L2, p=p,
                     design=protocol, Se=Se.H2, Sp=Sp.H2 )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           postGit=G, method="Bayes")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           emGit=G, emburn=burn, method="MAP")


p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 8)
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

## Mean:
Table.D6[i1+1:4, 1] <- round(colMeans(p.Mean), 3)   # Est
Table.D6[i1+1:4, 2] <- round(apply(p.Mean,2,sd), 4) # SD
Table.D6[i1+1:4, 3] <- round(colMeans(p.SE), 4)     # SE

## MAP:
Table.D6[i1+6:9, 1] <- round(colMeans(p.MAP), 3)    # Est
Table.D6[i1+6:9, 2] <- round(apply(p.MAP,2,sd), 4)  # SD
Table.D6[i1+6:9, 3] <- round(colMeans(p.SE), 4)     # SE

if( Configuration == "I" ){
  # Mean
  Table.D7[1:8, 1] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D7[1:8, 2] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D7[1:8, 3] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D7[10:17, 1] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D7[10:17, 2] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D7[10:17, 3] <- round(colMeans(delta.SE), 3)     # SE
}

if(Configuration == "II"){
  # Mean
  Table.D8[1:8, 1] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D8[1:8, 2] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D8[1:8, 3] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D8[10:17, 1] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D8[10:17, 2] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D8[10:17, 3] <- round(colMeans(delta.SE), 3)     # SE
}

## Avg. number of tests:
Table.D6[i1+11, 1] <- round( mean(T), 1 )    
Table.D6[i1+12, 1] <- round( sd(T), 1 )


#========= Estimation with H3 protocol =========#

protocol <- H3
list.zmat <- lapply( rep(N, nsims), hier.alg.data_L2, p=p,
                     design=protocol, Se=Se.H3, Sp=Sp.H3 )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           postGit=G, method="Bayes")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           emGit=G, emburn=burn, method="MAP")


p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 8)
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

## Mean:
Table.D6[i1+1:4, 5] <- round(colMeans(p.Mean), 3)   # Est
Table.D6[i1+1:4, 6] <- round(apply(p.Mean,2,sd), 4) # SD
Table.D6[i1+1:4, 7] <- round(colMeans(p.SE), 4)     # SE

## MAP:
Table.D6[i1+6:9, 5] <- round(colMeans(p.MAP), 3)    # Est
Table.D6[i1+6:9, 6] <- round(apply(p.MAP,2,sd), 4)  # SD
Table.D6[i1+6:9, 7] <- round(colMeans(p.SE), 4)     # SE

if( Configuration == "I" ){
  # Mean
  Table.D7[1:8, 5] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D7[1:8, 6] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D7[1:8, 7] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D7[10:17, 5] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D7[10:17, 6] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D7[10:17, 7] <- round(colMeans(delta.SE), 3)     # SE
}

if(Configuration == "II"){
  # Mean
  Table.D8[1:8, 5] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D8[1:8, 6] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D8[1:8, 7] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D8[10:17, 5] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D8[10:17, 6] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D8[10:17, 7] <- round(colMeans(delta.SE), 3)     # SE
}

## Avg. number of tests:
Table.D6[i1+11, 5] <- round( mean(T), 1 )    
Table.D6[i1+12, 5] <- round( sd(T), 1 )


#========= Estimation with H4 protocol =========#

protocol <- H4
list.zmat <- lapply( rep(N, nsims), hier.alg.data_L2, p=p,
                     design=protocol, Se=Se.H4, Sp=Sp.H4 )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           postGit=G, method="Bayes")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           emGit=G, emburn=burn, method="MAP")


p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 8)
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

## Mean:
Table.D6[i1+1:4, 9] <- round(colMeans(p.Mean), 3)   # Est
Table.D6[i1+1:4, 10] <- round(apply(p.Mean,2,sd), 4) # SD
Table.D6[i1+1:4, 11] <- round(colMeans(p.SE), 4)     # SE

## MAP:
Table.D6[i1+6:9, 9] <- round(colMeans(p.MAP), 3)    # Est
Table.D6[i1+6:9, 10] <- round(apply(p.MAP,2,sd), 4)  # SD
Table.D6[i1+6:9, 11] <- round(colMeans(p.SE), 4)     # SE

if( Configuration == "I" ){
  # Mean
  Table.D7[1:8, 9] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D7[1:8, 10] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D7[1:8, 11] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D7[10:17, 9] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D7[10:17, 10] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D7[10:17, 11] <- round(colMeans(delta.SE), 3)     # SE
}

if(Configuration == "II"){
  # Mean
  Table.D8[1:8, 9] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D8[1:8, 10] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D8[1:8, 11] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D8[10:17, 9] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D8[10:17, 10] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D8[10:17, 11] <- round(colMeans(delta.SE), 3)     # SE
}

## Avg. number of tests:
Table.D6[i1+11, 9] <- round( mean(T), 1 )    
Table.D6[i1+12, 9] <- round( sd(T), 1 )


#========= Estimation with AT protocol =========#

protocol <- AT
list.zmat <- lapply( rep(N, nsims), array.2dim.data_L2, p=p,
                     design=protocol, Se=Se.AT, Sp=Sp.AT )

## POSTERIOR SAMPLING 
res.bayes <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           postGit=G, method="Bayes")

## MAP ESTIMATION
res.map <- parLapply(cl, list.zmat, mult.gt.bayes_L2, 
           p0=p0, delta0=delta0, N=N, S=length(protocol), 
           p.pr=p.pr, se.pr=se.pr, sp.pr=sp.pr, 
           emGit=G, emburn=burn, method="MAP")


p.Mean <- p.MAP <- p.SE <- matrix(-9, nsims, 4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9, nsims, 8)
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

## Mean:
Table.D6[i1+1:4, 13] <- round(colMeans(p.Mean), 3)   # Est
Table.D6[i1+1:4, 14] <- round(apply(p.Mean,2,sd), 4) # SD
Table.D6[i1+1:4, 15] <- round(colMeans(p.SE), 4)     # SE

## MAP:
Table.D6[i1+6:9, 13] <- round(colMeans(p.MAP), 3)    # Est
Table.D6[i1+6:9, 14] <- round(apply(p.MAP,2,sd), 4)  # SD
Table.D6[i1+6:9, 15] <- round(colMeans(p.SE), 4)     # SE

if( Configuration == "I" ){
  # Mean
  Table.D7[1:8, 13] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D7[1:8, 14] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D7[1:8, 15] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D7[10:17, 13] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D7[10:17, 14] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D7[10:17, 15] <- round(colMeans(delta.SE), 3)     # SE
}

if(Configuration == "II"){
  # Mean
  Table.D8[1:8, 13] <- round(colMeans(delta.Mean), 3)   # Est
  Table.D8[1:8, 14] <- round(apply(delta.Mean,2,sd), 3) # SD
  Table.D8[1:8, 15] <- round(colMeans(delta.SE), 3)     # SE

  # MAP
  Table.D8[10:17, 13] <- round(colMeans(delta.MAP), 3)    # Est
  Table.D8[10:17, 14] <- round(apply(delta.MAP,2,sd), 3)  # SD
  Table.D8[10:17, 15] <- round(colMeans(delta.SE), 3)     # SE
}

## Avg. number of tests:
Table.D6[i1+11, 13] <- round( mean(T), 1 )    
Table.D6[i1+12, 13] <- round( sd(T), 1 )

}

stopCluster( cl )

## Preparing & saving Table D.6
dfr1 <- data.frame( cbind(
     c("","","","","","Config-I","","","","","","", 
	   "",  "","","","","","Config-II","","","","","",""),
     c("", "", "", "p00=0.95", "p10=0.02", "p01=0.02", "p11=0.01", "", "", "", "", "", "", 
	   "", "", "", "p00=0.990", "p10=0.004", "p01=0.004", "p11=0.002", "", "", "", "", ""),  
     c("Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T", "", 
	   "Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T")))

dfr1 <- data.frame( dfr1, Table.D6 )
dfr1 <- sapply(dfr1, as.character)
dfr1[is.na(dfr1)] <- " "
Table_D6 <- as.data.frame(dfr1)
colnames( Table_D6 ) <- c("", "", "", "Est", "SD", "SE", "", "Est", "SD", "SE", 
                           "", "Est", "SD", "SE", "", "Est", "SD", "SE")				   

saveRDS(Table_D6, "./simulation/intermediate_results/Table_D6.rds")


## Preparing & saving Table D.7 & D.8
dfr2 <- data.frame( cbind(
     c("", "", "", "", "Se:(1)1=0.95", "Se:(1)2=0.95", "Sp:(1)1=0.98", "Sp:(1)2=0.98",
       "Se:(2)1=0.98", "Se:(2)2=0.98", "Sp:(2)1=0.99", "Sp:(2)2=0.99", "", "", "", "", ""),
     c("", "", "", "Mean", "", "", "", "", "", "", "", "", "", "MAP", "", "", "")))

dfr3 <- data.frame( dfr2, Table.D7 )
dfr3 <- sapply(dfr3, as.character)
dfr3[is.na(dfr3)] <- " "

Table_D7 <- as.data.frame(dfr3)
colnames( Table_D7 ) <- c("", "", "Est", "SD", "SE", "", "Est", "SD", "SE", 
                           "", "Est", "SD", "SE", "", "Est", "SD", "SE")
						   
saveRDS(Table_D7, "./simulation/intermediate_results/Table_D7.rds")


dfr4 <- data.frame( dfr2, Table.D8 )
dfr4 <- sapply(dfr4, as.character)
dfr4[is.na(dfr4)] <- " "

Table_D8 <- as.data.frame(dfr4)
colnames( Table_D8 ) <- c("", "", "Est", "SD", "SE", "", "Est", "SD", "SE", 
                           "", "Est", "SD", "SE", "", "Est", "SD", "SE")
						   
saveRDS(Table_D8, "./simulation/intermediate_results/Table_D8.rds")

