
##################################################
# This document provides simulation evidence of 
# the article, "Estimating the prevalence of 
# two or more diseases using outcomes from multiplex 
# group testing," by Md S. Warasi, Joshua M. Tebbs,
# Christopher S. McMahan, and Christopher R. Bilder.
#
# The simulation uses: 
# (a) one assay, i.e., L=1;
# (b) KNOWN assay accuracy probabilities (Se & Sp);
# (c) flat Dirichlet prior on p=(p00,p10,p01,p11).
# 
# The simulation description is provided in Section 5.
#
# The results are shown in web-based Supporting Information.
# 
# Last updated on 12/14/2022.
##################################################

## Set the working directory:
setwd(dir = "C:/Users/msarker/Desktop/GitHub-MultiplexBayes_121422")

## Importing necessary functions:
source("MainFunction.txt")
source("SupportPrograms.txt")
source("MultStageHierData.txt")
source("TwoStageArrayData.txt")

## Package MCMCpack is required:
# install.packages("MCMCpack")

##################################################

## Set up simulation parameters here.

N <- 5000         # Sample size
Se <- c(.95,.95)  # True assay sensitivities
Sp <- c(.99,.99)  # True assay specificities

# Choose the disease prevalence configuration.
# We use configurations I & II in the paper.
# Configuration <- "I"   
Configuration <- "II"

if(Configuration == "I"){
  p <- c(.95,.02,.02,.01)     # True value of p
  IPT <- 5
  H2 <- c(5,1)
  H3 <- c(9,3,1)
  H4 <- c(18,6,3,1)
  AT <- c(11,11,1)  # For AT, c(11,11,1), not c(11,11)
}                   # because tests are done in 3 stages:
                    # row pool, col pool, & individual testing

if(Configuration == "II"){
  p <- c(.990,.004,.004,.002) # True value of p
  IPT <- 11
  H2 <- c(11,1)
  H3 <- c(25,5,1)
  H4 <- c(48,12,4,1)
  AT <- c(29,29,1)
}

## Specify the prior hyperparameters. 
## We use flat priors on all parameters as shown below.
p.pr <- rep(1,4)

## Choose the number of simulated data sets.
## We use sims = 500 data sets in the paper.
## Note: Computing time with 500 could be long.
sims <- 500

G <- 12000               # number of Gibbs iterates
burn <- 2000             # a burn-in period of 2000
pick <- seq(1,10000,5)   # keeping every 5th

## Choose an initial value. One can start at the true 
## value, the default value p0=c(.90,.06,.03,.01), 
## or any reasonable choice of p: 
p0 <- p                  


## Computing starts here under each testing protocol >>


#========= Simulation using IPT protocol =========#

p.Mean <- p.MAP <- p.SE <- matrix(-9,sims,4)
T <- rep(-9,sims)

for(s in 1:sims){

## Simulating data
protocol <- IPT
out <- hier.alg.data(p=p,N=N,design=protocol,Se=Se,Sp=Sp)
Z <- out$Data    
T[s] <- out$T

## POSTERIOR SAMPLING ALGORITHM    
res1 <- mult.gt.bayes(p0=p0,
           Z=Z,N=N,S=length(protocol),p.pr=p.pr,
           postGit=G,method="Bayes",accuracy="known")
p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION VIA EM ALGORITHM
res2 <- mult.gt.bayes(p0=p0,Z=Z,N=N,S=length(protocol),
                p.pr=p.pr,emGit=G,emburn=burn,method="MAP",
                accuracy="known")
p.MAP[s, ] <- res2$prevalence

## Printing the number of data sets completed:
print(s)

}

## Mean Estimation:
round(colMeans(p.Mean),3)   # Est
round(apply(p.Mean,2,sd),4) # SD
round(colMeans(p.SE),4)     # SE

## MAP Estimation:
round(colMeans(p.MAP),3)    # Est
round(apply(p.MAP,2,sd),4)  # SD
round(colMeans(p.SE),4)     # SE

## Avg. number of tests:
round( mean(T), 1 )    



#========= Simulation using H2 protocol =========#

p.Mean <- p.MAP <- p.SE <- matrix(-9,sims,4)
T <- rep(-9,sims)

for(s in 1:sims){

## Simulating data
protocol <- H2
out <- hier.alg.data(p=p,N=N,design=protocol,Se=Se,Sp=Sp)
Z <- out$Data    
T[s] <- out$T

## POSTERIOR SAMPLING ALGORITHM    
res1 <- mult.gt.bayes(p0=p0,
           Z=Z,N=N,S=length(protocol),p.pr=p.pr,
           postGit=G,method="Bayes",accuracy="known")
p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION VIA EM ALGORITHM
res2 <- mult.gt.bayes(p0=p0,Z=Z,N=N,S=length(protocol),
                p.pr=p.pr,emGit=G,emburn=burn,method="MAP",
                accuracy="known")
p.MAP[s, ] <- res2$prevalence

## Printing the number of data sets completed:
print(s)

}

## Mean Estimation:
round(colMeans(p.Mean),3)   # Est
round(apply(p.Mean,2,sd),4) # SD
round(colMeans(p.SE),4)     # SE

## MAP Estimation:
round(colMeans(p.MAP),3)    # Est
round(apply(p.MAP,2,sd),4)  # SD
round(colMeans(p.SE),4)     # SE

## Avg. number of tests:
round( mean(T), 1 )    


#========= Simulation using H3 protocol =========#

p.Mean <- p.MAP <- p.SE <- matrix(-9,sims,4)
T <- rep(-9,sims)

for(s in 1:sims){

## Simulating data
protocol <- H3
out <- hier.alg.data(p=p,N=N,design=protocol,Se=Se,Sp=Sp)
Z <- out$Data    
T[s] <- out$T

## POSTERIOR SAMPLING ALGORITHM    
res1 <- mult.gt.bayes(p0=p0,
           Z=Z,N=N,S=length(protocol),p.pr=p.pr,
           postGit=G,method="Bayes",accuracy="known")
p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION VIA EM ALGORITHM
res2 <- mult.gt.bayes(p0=p0,Z=Z,N=N,S=length(protocol),
                p.pr=p.pr,emGit=G,emburn=burn,method="MAP",
                accuracy="known")
p.MAP[s, ] <- res2$prevalence

## Printing the number of data sets completed:
print(s)

}

## Mean Estimation:
round(colMeans(p.Mean),3)   # Est
round(apply(p.Mean,2,sd),4) # SD
round(colMeans(p.SE),4)     # SE

## MAP Estimation:
round(colMeans(p.MAP),3)    # Est
round(apply(p.MAP,2,sd),4)  # SD
round(colMeans(p.SE),4)     # SE

## Avg. number of tests:
round( mean(T), 1 )    


#========= Simulation using H4 protocol =========#

p.Mean <- p.MAP <- p.SE <- matrix(-9,sims,4)
T <- rep(-9,sims)

for(s in 1:sims){

## Simulating data
protocol <- H4
out <- hier.alg.data(p=p,N=N,design=protocol,Se=Se,Sp=Sp)
Z <- out$Data    
T[s] <- out$T

## POSTERIOR SAMPLING ALGORITHM    
res1 <- mult.gt.bayes(p0=p0,
           Z=Z,N=N,S=length(protocol),p.pr=p.pr,
           postGit=G,method="Bayes",accuracy="known")
p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION VIA EM ALGORITHM
res2 <- mult.gt.bayes(p0=p0,Z=Z,N=N,S=length(protocol),
                p.pr=p.pr,emGit=G,emburn=burn,method="MAP",
                accuracy="known")
p.MAP[s, ] <- res2$prevalence

## Printing the number of data sets completed:
print(s)

}

## Mean Estimation:
round(colMeans(p.Mean),3)   # Est
round(apply(p.Mean,2,sd),4) # SD
round(colMeans(p.SE),4)     # SE

## MAP Estimation:
round(colMeans(p.MAP),3)    # Est
round(apply(p.MAP,2,sd),4)  # SD
round(colMeans(p.SE),4)     # SE

## Avg. number of tests:
round( mean(T), 1 )    


#========= Simulation using AT protocol =========#

p.Mean <- p.MAP <- p.SE <- matrix(-9,sims,4)
T <- rep(-9,sims)

for(s in 1:sims){

## Simulating data
protocol <- AT
out <- array.2dim.data(p=p,N=N,design=protocol,Se=Se,Sp=Sp)
Z <- out$Data    
T[s] <- out$T

## POSTERIOR SAMPLING ALGORITHM    
res1 <- mult.gt.bayes(p0=p0,
           Z=Z,N=N,S=length(protocol),p.pr=p.pr,
           postGit=G,method="Bayes",accuracy="known")
p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION VIA EM ALGORITHM
res2 <- mult.gt.bayes(p0=p0,Z=Z,N=N,S=length(protocol),
                p.pr=p.pr,emGit=G,emburn=burn,method="MAP",
                accuracy="known")
p.MAP[s, ] <- res2$prevalence

## Printing the number of data sets completed:
print(s)

}

## Mean Estimation:
round(colMeans(p.Mean),3)   # Est
round(apply(p.Mean,2,sd),4) # SD
round(colMeans(p.SE),4)     # SE

## MAP Estimation:
round(colMeans(p.MAP),3)    # Est
round(apply(p.MAP,2,sd),4)  # SD
round(colMeans(p.SE),4)     # SE

## Avg. number of tests:
round( mean(T), 1 )    

