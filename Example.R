
##############################################

## In this R script, we provide an example, 
## which is used to the README.md file.

##############################################

## Set the working directory:
setwd(dir = "C:\programs")

# Load the package and code:
library(MCMCpack)
source("MainFunctions.R")
source("SupportPrograms.R")

## Reading in the data:
Z <- read.csv( "Urine_GT_H3.csv", row.names=1 )
Z <- as.matrix(Z)

## Sample size:
N <- 4402

## Three-stage hierarchical protocol with pool sizes:
protocol <- c(6, 2, 1)

head( Z )

tail(Z)


## Specify the prior hyperparameters as shown below. 
## Flat priors have been used for all parameters.
p.pr <- rep(1, 4)
Se1.pr <- c(1, 1)
Se2.pr <- c(1, 1)
Sp1.pr <- c(1, 1)
Sp2.pr <- c(1, 1)

G <- 12000                # number of Gibbs iterates
burn <- 2000              # burn-in period
pick <- seq(1, 10000, 5)  # keeping every 5th

## Convergence tolerance for the EM algorithm
epsilon <- 0.001

## Maximum number of iterations for the EM algorithm
emmaxit <- 200

## Specify an initial value for p and delta. These values can be 
## specified from historical, pilot study, or any contextual estimate. 
## We use p0=c(.92,.05,.02,.01) and delta0=c(.96,.96,.98,.98) 
## throughout the article. 
p0 <- c(.92, .05, .02, .01)
delta0 <- c(.96, .96, .98, .98)

set.seed( 123 )


## Case I: KNOWN accuracy probabilities

## POSTERIOR SAMPLING 
post.samp <- mult.gt.bayes( Z=Z, p0=p0, N=N, S=length(protocol), 
                            p.pr=p.pr, postGit=G, method="Bayes", 
                            accuracy="known")

# Mean estimate:
round( colMeans(post.samp$prevalence[-(1:burn), ][pick, ]), 4)
round( apply(post.samp$prevalence[-(1:burn), ][pick, ], 2, sd), 4)


## MAP ESTIMATION
res.map <- mult.gt.bayes( Z=Z, p0=p0, N=N, S=length(protocol),
                          p.pr=p.pr, emGit=G, emburn=burn, 
                          method="MAP", accuracy="known" )

round( res.map$prevalence, 4)

## Number of tests expended:
nrow( Z )



## Case II: UNKNOWN accuracy probabilities

Z[ ,4:7] <- -9

head( Z )

tail( Z )


## POSTERIOR SAMPLING 
post.samp <- mult.gt.bayes( Z=Z, p0=p0, delta0=delta0, N=N, 
                 S=length(protocol), p.pr=p.pr, Se1.pr=Se1.pr, Se2.pr=Se2.pr,
                 Sp1.pr=Sp1.pr, Sp2.pr=Sp2.pr, postGit=G, method="Bayes", 
                 accuracy="unknown" )

# Mean estimate:
round( colMeans(post.samp$prevalence[-(1:burn), ][pick, ]), 4)
round( apply(post.samp$prevalence[-(1:burn), ][pick, ], 2, sd), 4)

round( colMeans(post.samp$accuracy[-(1:burn), ][pick, ]), 4)
round( apply(post.samp$accuracy[-(1:burn), ][pick, ], 2, sd), 4)


## MAP ESTIMATION
res.map <-  mult.gt.bayes( Z=Z, p0=p0, delta0=delta0, N=N, 
                 S=length(protocol), p.pr=p.pr, Se1.pr=Se1.pr, Se2.pr=Se2.pr, 
                 Sp1.pr=Sp1.pr, Sp2.pr=Sp2.pr, emGit=G, emburn=burn, 
                 method="MAP", accuracy="unknown" )
			
round( res.map$prevalence, 4)

round( res.map$accuracy, 4)

## Number of tests expended:
nrow( Z )


			
			