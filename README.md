# General-MultiplexBayes

This repository contains R programs of the article, "Estimating the prevalence of two or more diseases using outcomes from multiplex group testing." Two main R functions (mult.gt.bayes for L=1 assay and mult.gt.bayes_L2 for L=2 assays) are provided, each of which can implement the posterior sampling algorithm and the EM algorithm proposed in Warasi et al. (2022+). This article uses multiplex group testing data for estimating the coinfection probabilities and the assay accuracy probabilities (sensitivity and specificity).

R code of the simulation examples in the article is split into three files:

1. Simulation1.R -- to reproduce the estimation results shown in Table D.1 (in Web Appendix D).

2. Simulation2.R -- to reproduce the estimation results shown in Tables 2-3 (in the article).

3. Simulation3.R -- to reproduce the estimation results shown in Tables D.3-D.8 (in Web Appendix D). 



Reference

Warasi, M., Tebbs, J., McMahan, C., and Bilder, C. (2023+). Estimating the prevalence of two or more diseases using outcomes from multiplex group testing. Under review.


##################################################

We herein show how the estimation techniques in Warasi et al. (2023+) can be implemented using the R functions provided in this repository. 
We illustrate the R functions using simulated data. The data are generated using the parameter configuration described in Section 5 (Tables 1-3).

##################################################

## Set the working directory:
setwd(dir = "C:\\programs")

## Importing necessary functions:
source("MainFunction.txt")

source("SupportPrograms.txt")

source("MultStageHierData.txt")

source("TwoStageArrayData.txt")

## Package MCMCpack is required:

install.packages("MCMCpack")

##################################################

### Data simulation using function: hier.alg.data()

### Parameter configurations for data simulation:
N <- 5000                       # Sample size

H3 <- c(9, 3, 1)                # H3 protocol

p <- c(0.95, 0.02, 0.02, 0.01)  # True parameter (Configuration I)

Se <- c(.95,.95)                # True Se for diseases 1 & 2

Sp <- c(.99,.99)                # True Sp for diseases 1 & 2

## Simulating data 
set.seed( 123 )

protocol <- H3

out <- hier.alg.data(p=p,N=N,design=protocol,Se=Se,Sp=Sp)

Z <- out$Data        # Simulated data is Z

T <- out$T           # The number of tests expended


The simulated data are in the matrix object Z. The data with other information need to be structured in a specific manner. The first two columns consist of test responses for diseases 1 & 2, the third column consists of pool sizes, columns 4-5 must have the assay sensitivities for diseases 1 & 2, columns 6-7 consist of specificities for diseases 1 & 2, and column 8 onward must have identification numbers of the individuals assigned to each pool. The pool identification numbers (i.e., row names) are not used so they do not need to be specified. Parts of a simulated data set is show in "MainFunction.txt" in this repository. For more information, please refer to the simulation examples provided in Simulation1.R, Simulation2.R, and Simulation3.R.

head( Z )

tail( Z )


##################################################

### Model fitting using function: mult.gt.bayes(). Specify the prior hyperparameters. We use flat priors on all parameters as shown below.

p.pr <- rep(1,4)

Se1.pr <- c(1,1)

Se2.pr <- c(1,1)

Sp1.pr <- c(1,1)

Sp2.pr <- c(1,1)

G <- 12000               # number of Gibbs iterates

burn <- 2000             # a burn-in period of 2000

pick <- seq(1,10000,5)   # keeping every 5th

## Choose an initial value. One can start at the true value, the default value p0=c(.90,.06,.03,.01), delta0=c(.95,.95,.98,.98), or any reasonable choice: 

p0 <- c(.90,.06,.03,.01)    

delta0 <- c(Se, Sp)

## POSTERIOR SAMPLING ALGORITHM    
res1 <- mult.gt.bayes(p0=p0,delta0=delta0,
           Z=Z,N=N,S=length(protocol),p.pr=p.pr,
           Se1.pr=Se1.pr,Se2.pr=Se2.pr,Sp1.pr=Sp1.pr,Sp2.pr=Sp2.pr,
           postGit=G,method="Bayes",accuracy="unknown")

## Mean estimation results are shown below. This computing is completed in 21 seconds on a computer that has an Intel 2.60 GHz processor and 32 GB of RAM.

p.Mean <- colMeans(res1$prevalence[-(1:burn),][pick,])

> p.Mean
   
0.9518  0.0191  0.0193  0.0098  

p.SE <- apply(res1$prevalence[-(1:burn),][pick,],2,sd)

> p.SE
   
0.0034  0.0022  0.0022  0.0014  

delta.Mean <- colMeans(res1$accuracy[-(1:burn),][pick,])

> delta.Mean

0.9365   0.9499   0.9943   0.9914

delta.SE <- apply(res1$accuracy[-(1:burn),][pick,],2,sd)

> delta.SE

0.0169   0.0140   0.0030   0.0033


## MAP ESTIMATION VIA EM ALGORITHM
res2 <- mult.gt.bayes(p0=p0,delta0=delta0,
            Z=Z,N=N,S=length(protocol),p.pr=p.pr,
            Se1.pr=Se1.pr,Se2.pr=Se2.pr,Sp1.pr=Sp1.pr,Sp2.pr=Sp2.pr,
            emGit=G,emburn=burn,method="MAP",accuracy="unknown")

## MAP estimation results are shown below. This computing is completed in 51 seconds that has an Intel 2.60 GHz processor and 32 GB of RAM.

p.MAP <- res2$prevalence

> p.MAP

0.9522  0.0190  0.0191  0.0097

delta.MAP <- res2$accuracy

> delta.MAP

0.9371  0.9521  0.9957  0.9925

