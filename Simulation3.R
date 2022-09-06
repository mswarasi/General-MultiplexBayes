
##################################################
# This document provides simulation evidence of 
# the article, "Estimating the prevalence of 
# two or more diseases using outcomes from multiplex 
# group testing," by Md S. Warasi, Joshua M. Tebbs,
# Christopher S. McMahan, and Christopher R. Bilder.
#
# The simulation uses: 
# (a) two assays, i.e., L=2;
# (b) UNKNOWN assay accuracy probabilities (Se & Sp);
# (c) flat Dirichlet prior on p=(p00,p10,p01,p11)
#     and but both flat and informative beta priors 
#     on Se's and Sp's.
# 
# The simulation description is provided in Section 5
# and web-based Supporting Information.
#
# The results are shown in web-based Supporting Information.
# 
# Last updated on 08/28/2022.
##################################################

rm(list=ls())

## Set the working directory:
setwd(dir = "C:\\programs")

## Importing necessary functions:
source("MainFunction_L2.txt")
source("SupportPrograms_L2.txt")
source("MultStageHierData_L2.txt")
source("TwoStageArrayData_L2.txt")

## Package MCMCpack is required:
# install.packages("MCMCpack")

##################################################

## Set up simulation parameters here.

N <- 5000         # Sample size

# For H2 simulation
Se.H2 <- rbind(c(.95,.95), 
			   c(.98,.98) )
Sp.H2 <- rbind(c(.98,.98), 
			   c(.99,.99) )
			  
# For H3 simulation
Se.H3 <- rbind(c(.95,.95), 
               c(.95,.95), 
			   c(.98,.98) )
Sp.H3 <- rbind(c(.98,.98), 
               c(.98,.98), 
			   c(.99,.99) )

# For H4 simulation
Se.H4 <- rbind(c(.95,.95), 
               c(.95,.95), 
               c(.95,.95),
			   c(.98,.98) )
Sp.H4 <- rbind(c(.98,.98), 
               c(.98,.98), 
               c(.98,.98), 
		  	   c(.99,.99) )

# For AT simulation
Se.AT <- rbind(c(.95,.95), 
			   c(.98,.98) )
Sp.AT <- rbind(c(.98,.98), 
			   c(.99,.99) )

## Choose the disease prevalence configuration.
# We use configurations I & II in the paper.
Configuration <- "I"   
# Configuration <- "II"

if(Configuration == "I"){
  p <- c(.95,.02,.02,.01)  # True value of p
  H2 <- c(5,1)
  H3 <- c(9,3,1)
  H4 <- c(18,6,3,1)
  AT <- c(11,11,1)  # For AT, c(11,11,1), not c(11,11)
}                   # because tests are done in 3 stages:
                    # row pool, col pool, & individual testing

if(Configuration == "II"){
  p <- c(.990,.004,.004,.002) # True value of p
  H2 <- c(11,1)
  H3 <- c(24,6,1)
  H4 <- c(36,12,4,1)
  AT <- c(29,29,1)
}

## True value of delta that we are estimating:
delta <- c(rep(.95,2),rep(.98,2),rep(.98,2),rep(.99,2))


## Specify the prior hyperparameters.
# Flat prior on p.
p.pr <- rep(1,4)

# Either flat or informative priors on Se & Sp:
SeSp.pr <- "flat" 
#SeSp.pr <- "informative" 

if( SeSp.pr == "flat" ){
  se.pr <- matrix(1,2,4)
  sp.pr <- matrix(1,2,4)
}
if( SeSp.pr == "informative" ){
  se.pr <- rbind(c(109.0, 6.7, 109.0, 6.7),  # for stage 1
                 c(100.0, 3.0, 100.0, 3.0))  # for stage 2

  sp.pr <- rbind(c(100.0, 3.0, 100.0, 3.0),  # for stage 1
                 c( 55.2, 1.6,  55.2, 1.6))  # for stage 2
}

## Choose the number of simulated data sets.
## We use sims = 500 data sets in the paper.
## Note: Computing time with 500 could be long.
sims <- 500

G <- 12000               # number of Gibbs iterates
burn <- 2000             # a burn-in period of 2000
pick <- seq(1,10000,5)   # keeping every 5th

## Choose an initial value. One can start at the true 
## value, the default value, or any reasonable choice: 
p0 <- p
delta0 <- c(rep(.95,2),rep(.98,2),rep(.98,2),rep(.99,2))

## Computing starts here under each testing protocol >>


#========= Simulation using H2 protocol =========#

p.Mean <- p.MAP <- p.SE <- matrix(-9,sims,4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9,sims,8)
T <- rep(-9,sims)

for(s in 1:sims){

## Simulating data
protocol <- H2
out <- hier.alg.data_L2(p=p,N=N,design=protocol,Se=Se.H2,Sp=Sp.H2)
zmat <- out$Data
T[s] <- out$T

## POSTERIOR SAMPLING ALGORITHM    
res1 <- mult.gt.bayes_L2(p0=p0,delta0=delta0,Z=zmat,
           N=N,S=length(protocol),p.pr=p.pr,se.pr=se.pr,
           postGit=G,method="Bayes")

p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)
delta.Mean[s, ] <- colMeans(res1$accuracy[-(1:burn),][pick,])
delta.SE[s, ] <- apply(res1$accuracy[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION VIA EM ALGORITHM
res2 <- mult.gt.bayes_L2(p0=p0,delta0=delta0,Z=zmat,
           N=N,S=length(protocol),p.pr=p.pr,se.pr=se.pr,
           sp.pr=sp.pr,emGit=G,emburn=burn,method="MAP")
p.MAP[s, ] <- res2$prevalence
delta.MAP[s, ] <- res2$accuracy

## Printing the number of data sets completed:
print(s)

}

## True parameters: 
p
delta


## Mean Estimation:
round(colMeans(p.Mean),3)   # Est
round(apply(p.Mean,2,sd),4) # SD
round(colMeans(p.SE),4)     # SE

round(colMeans(delta.Mean),3)   # Est
round(apply(delta.Mean,2,sd),4) # SD
round(colMeans(delta.SE),4)     # SE

## MAP Estimation:
round(colMeans(p.MAP),3)    # Est
round(apply(p.MAP,2,sd),4)  # SD
round(colMeans(p.SE),4)     # SE

round(colMeans(delta.MAP),3)    # Est
round(apply(delta.MAP,2,sd),4)  # SD
round(colMeans(delta.SE),4)     # SE

## Avg. number of tests:
round( mean(T), 1 )    


#========= Simulation using H3 protocol =========#

p.Mean <- p.MAP <- p.SE <- matrix(-9,sims,4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9,sims,8)
T <- rep(-9,sims)

for(s in 1:sims){

## Simulating data
protocol <- H3
out <- hier.alg.data_L2(p=p,N=N,design=protocol,Se=Se.H3,Sp=Sp.H3)
zmat <- out$Data
T[s] <- out$T

## POSTERIOR SAMPLING ALGORITHM    
res1 <- mult.gt.bayes_L2(p0=p0,delta0=delta0,Z=zmat,
           N=N,S=length(protocol),p.pr=p.pr,se.pr=se.pr,
           postGit=G,method="Bayes")

p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)
delta.Mean[s, ] <- colMeans(res1$accuracy[-(1:burn),][pick,])
delta.SE[s, ] <- apply(res1$accuracy[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION VIA EM ALGORITHM
res2 <- mult.gt.bayes_L2(p0=p0,delta0=delta0,Z=zmat,
           N=N,S=length(protocol),p.pr=p.pr,se.pr=se.pr,
           sp.pr=sp.pr,emGit=G,emburn=burn,method="MAP")
p.MAP[s, ] <- res2$prevalence
delta.MAP[s, ] <- res2$accuracy

## Printing the number of data sets completed:
print(s)

}

## True parameters: 
p
delta


## Mean Estimation:
round(colMeans(p.Mean),3)   # Est
round(apply(p.Mean,2,sd),4) # SD
round(colMeans(p.SE),4)     # SE

round(colMeans(delta.Mean),3)   # Est
round(apply(delta.Mean,2,sd),4) # SD
round(colMeans(delta.SE),4)     # SE

## MAP Estimation:
round(colMeans(p.MAP),3)    # Est
round(apply(p.MAP,2,sd),4)  # SD
round(colMeans(p.SE),4)     # SE

round(colMeans(delta.MAP),3)    # Est
round(apply(delta.MAP,2,sd),4)  # SD
round(colMeans(delta.SE),4)     # SE

## Avg. number of tests:
round( mean(T), 1 )


#========= Simulation using H4 protocol =========#

p.Mean <- p.MAP <- p.SE <- matrix(-9,sims,4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9,sims,8)
T <- rep(-9,sims)

for(s in 1:sims){

## Simulating data
protocol <- H4
out <- hier.alg.data_L2(p=p,N=N,design=protocol,Se=Se.H4,Sp=Sp.H4)
zmat <- out$Data
T[s] <- out$T

## POSTERIOR SAMPLING ALGORITHM    
res1 <- mult.gt.bayes_L2(p0=p0,delta0=delta0,Z=zmat,
           N=N,S=length(protocol),p.pr=p.pr,se.pr=se.pr,
           postGit=G,method="Bayes")

p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)
delta.Mean[s, ] <- colMeans(res1$accuracy[-(1:burn),][pick,])
delta.SE[s, ] <- apply(res1$accuracy[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION VIA EM ALGORITHM
res2 <- mult.gt.bayes_L2(p0=p0,delta0=delta0,Z=zmat,
           N=N,S=length(protocol),p.pr=p.pr,se.pr=se.pr,
           sp.pr=sp.pr,emGit=G,emburn=burn,method="MAP")
p.MAP[s, ] <- res2$prevalence
delta.MAP[s, ] <- res2$accuracy

## Printing the number of data sets completed:
print(s)

}

## True parameters: 
p
delta


## Mean Estimation:
round(colMeans(p.Mean),3)   # Est
round(apply(p.Mean,2,sd),4) # SD
round(colMeans(p.SE),4)     # SE

round(colMeans(delta.Mean),3)   # Est
round(apply(delta.Mean,2,sd),4) # SD
round(colMeans(delta.SE),4)     # SE

## MAP Estimation:
round(colMeans(p.MAP),3)    # Est
round(apply(p.MAP,2,sd),4)  # SD
round(colMeans(p.SE),4)     # SE

round(colMeans(delta.MAP),3)    # Est
round(apply(delta.MAP,2,sd),4)  # SD
round(colMeans(delta.SE),4)     # SE

## Avg. number of tests:
round( mean(T), 1 )    


#========= Simulation using AT protocol =========#

p.Mean <- p.MAP <- p.SE <- matrix(-9,sims,4)
delta.Mean <- delta.MAP <- delta.SE <- matrix(-9,sims,8)
T <- rep(-9,sims)

for(s in 1:sims){

## Simulating data
protocol <- AT
out <- array.2dim.data_L2(p=p,N=N,design=protocol,Se=Se.AT,Sp=Sp.AT)
zmat <- out$Data
T[s] <- out$T

## POSTERIOR SAMPLING ALGORITHM    
res1 <- mult.gt.bayes_L2(p0=p0,delta0=delta0,Z=zmat,
           N=N,S=length(protocol),p.pr=p.pr,se.pr=se.pr,
           postGit=G,method="Bayes")

p.Mean[s, ] <- colMeans(res1$prevalence[-(1:burn),][pick,])
p.SE[s, ] <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)
delta.Mean[s, ] <- colMeans(res1$accuracy[-(1:burn),][pick,])
delta.SE[s, ] <- apply(res1$accuracy[-(1:burn),][pick,],2,sd)

## MAP ESTIMATION VIA EM ALGORITHM
res2 <- mult.gt.bayes_L2(p0=p0,delta0=delta0,Z=zmat,
           N=N,S=length(protocol),p.pr=p.pr,se.pr=se.pr,
           sp.pr=sp.pr,emGit=G,emburn=burn,method="MAP")
p.MAP[s, ] <- res2$prevalence
delta.MAP[s, ] <- res2$accuracy

## Printing the number of data sets completed:
print(s)

}

## True parameters: 
p
delta


## Mean Estimation:
round(colMeans(p.Mean),3)   # Est
round(apply(p.Mean,2,sd),4) # SD
round(colMeans(p.SE),4)     # SE

round(colMeans(delta.Mean),3)   # Est
round(apply(delta.Mean,2,sd),4) # SD
round(colMeans(delta.SE),4)     # SE

## MAP Estimation:
round(colMeans(p.MAP),3)    # Est
round(apply(p.MAP,2,sd),4)  # SD
round(colMeans(p.SE),4)     # SE

round(colMeans(delta.MAP),3)    # Est
round(apply(delta.MAP,2,sd),4)  # SD
round(colMeans(delta.SE),4)     # SE

## Avg. number of tests:
round( mean(T), 1 )    

