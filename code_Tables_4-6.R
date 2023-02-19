
##################################################

# 'ncores' is the number of cores to be used
# for parallel computing in Windows with 64-bit R:
ncores <- 10

# 'nsims' is the number of data replicates:
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
## shown in Table 4.

## The sensitivities and specificities are obtained
## from the Aptima Combo 2 Assay package inserts.
S <- 3
se.urine <- cbind( rep(0.947, S), rep(0.913, S) )
sp.urine <- cbind( rep(0.989, S), rep(0.993, S) )

se.swab <- cbind( rep(0.942, S), rep(0.992, S) )
sp.swab <- cbind( rep(0.976, S), rep(0.987, S) )

## The vectors of coinfection probabilities shown below  
## are obtained from the individual case identification 
## outcomes in the Iowa CT/NG data. These are viewed as 
## historical estimates of p when calculating the optimal 
## pool sizes presented in Table 4 (Section 6).
p.urine <- c( 0.9082, 0.0811, 0.0057, 0.0050 )
p.swab <-  c( 0.9086, 0.0812, 0.0054, 0.0048 )

opt.psz <- array(NA, c(S, S, 2), )
stratum <- c( "urine", "swab" )

for( d in 1:2 ){

  if( stratum[d] == "urine" ){ ## urine stratum
    p <- p.urine
    se.mat <- se.urine
    sp.mat <- sp.urine
  }
  if( stratum[d] == "swab" ){  ## swab stratum
    p <- p.swab
    se.mat <- se.swab
    sp.mat <- sp.swab
  }

  ## Hierarchical testing
  for(s in 2:S){
    designs <- read.table( paste("./source_code/","psz_s", s, ".csv", sep="") )
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
  max.ary.sz <- 20    ## max array size
  designs <- 2:max.ary.sz
  at.et <- rep(-9, length(designs))
  for(k in designs){
    res <- OTC2(algorithm = "A2", p.vec = p, Se = Se[1:2, ], 
                Sp = Sp[1:2, ], group.sz = k, trace = FALSE, print.time = FALSE)
    at.et[k-1] <- res[5]$opt.ET$value
  }
  opt.psz[S, 1:2, d] <- designs[at.et == min(at.et)]
}

Table.4 <- cbind(
  c("Urine", "", "", "", "Swab", "", ""),
  c("H2", "H3", "AT", "", "H2", "H3", "AT"),
  c(opt.psz[1,1,1], opt.psz[2,1,1], opt.psz[3,1,1], "", 
    opt.psz[1,1,2], opt.psz[2,1,2], opt.psz[3,1,2]),	
  c(opt.psz[1,2,1], opt.psz[2,2,1], opt.psz[3,2,1], "", 
    opt.psz[1,2,2], opt.psz[2,2,2], opt.psz[3,2,2]),
  c(opt.psz[1,3,1], opt.psz[2,3,1], opt.psz[3,3,1], "", 
    opt.psz[1,3,2], opt.psz[2,3,2], opt.psz[3,3,2]),	
  c("", "", "", "", "", "", ""),
#  c("Urine", "", "", "", "Swab", "", ""),
  c("CT", "NG", "", "", "CT", "NG", ""),
  c( paste("Se:(1)1=", se.urine[1, 1], sep=""),  
  paste("Se:(1)2=", se.urine[1, 2], sep=""), "",
  "", paste("Se:(1)1=", se.swab[1, 1], sep=""),  
  paste("Se:(1)2=", se.swab[1, 2], sep=""), "" ),
  c( paste("Sp:(1)1=", sp.urine[1, 1], sep=""),  
  paste("Sp:(1)2=", sp.urine[1, 2], sep=""), "",
  "", paste("Sp:(1)1=", sp.swab[1, 1], sep=""),  
  paste("Sp:(1)2=", sp.swab[1, 2], sep=""), "" )
)

Table.4[ is.na(Table.4) ] <- ""
Table_4 <- as.data.frame( Table.4 )
colnames( Table_4 ) <- c( "", "Protocol", "Pool sizes", "", "", "", "", "Sensitivity", "Specificity" )
saveRDS(Table_4, "./data_analysis/intermediate_results/Table_4.rds")

##################################################

## The code below produces the results shown in Tables 5 and 6. 
## The intermediate results shown here are based on simulated
## data because the real data cannot be shared. 

data_urine <- readRDS("./data_analysis/Iowa_data/data_urine_stratum.rds")
data_swab <- readRDS("./data_analysis/Iowa_data/data_swab_stratum.rds")

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

## Remarks about initial values:
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

# Setting seed:
RNGkind(kind = "L'Ecuyer-CMRG")
set.seed( 123 )
clusterSetRNGStream(cl = cl, iseed = 123)

options( scipen = 999 )
Table.5 <- matrix(NA, nrow = 25, ncol = 8)
Table.6 <- matrix(NA, nrow = 19, ncol = 8)

for( stratum in c( "urine", "swab" ) ){

if( stratum == "urine" ){
  Se <- se.urine
  Sp <- sp.urine
  Yt <- data_urine
  N <- 4402
  H2 <- c(4, 1)
  H3 <- c(9, 3, 1)
  AT <- c(8, 8, 1)
  i1 <- 0
  i2 <- 0
}

if( stratum == "swab" ){
  Se <- se.swab
  Sp <- sp.swab
  Yt <- data_swab
  N <- 10048
  H2 <- c(4, 1)
  H3 <- c(9, 3, 1)
  AT <- c(8, 8, 1)
  i1 <- 13
  i2 <- 10
}

#========= Estimation with H2 protocol =========#

protocol <- H2
list.zmat <- lapply( rep(N, nsims), hier.alg.data, 
                     design=protocol, Se=Se, Sp=Sp, Yt=Yt )

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
Table.5[i1+1:4, 1] <- round(colMeans(p.Mean), 3)   # Est
Table.5[i1+1:4, 2] <- round(colMeans(p.SE), 4)     # SE

Table.6[i2+1:4, 1] <- round(colMeans(delta.Mean), 3)   # Est
Table.6[i2+1:4, 2] <- round(colMeans(delta.SE), 3)     # SE

## MAP Estimation:
Table.5[i1+6:9, 1] <- round(colMeans(p.MAP), 3)    # Est
Table.5[i1+6:9, 2] <- round(colMeans(p.SE), 4)     # SE

Table.6[i2+6:9, 1] <- round(colMeans(delta.MAP), 3)    # Est
Table.6[i2+6:9, 2] <- round(colMeans(delta.SE), 3)     # SE

## Avg. number of tests:
Table.5[i1+11, 1] <- round( mean(T), 1 )    
Table.5[i1+12, 1] <- round( sd(T), 1 )

#========= Estimation with H3 protocol =========#

protocol <- H3
list.zmat <- lapply( rep(N, nsims), hier.alg.data,
                     design=protocol, Se=Se, Sp=Sp, Yt=Yt )

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
Table.5[i1+1:4, 4] <- round(colMeans(p.Mean), 3)   # Est
Table.5[i1+1:4, 5] <- round(colMeans(p.SE), 4)     # SE

Table.6[i2+1:4, 4] <- round(colMeans(delta.Mean), 3)   # Est
Table.6[i2+1:4, 5] <- round(colMeans(delta.SE), 3)     # SE

## MAP Estimation:
Table.5[i1+6:9, 4] <- round(colMeans(p.MAP), 3)    # Est
Table.5[i1+6:9, 5] <- round(colMeans(p.SE), 4)     # SE

Table.6[i2+6:9, 4] <- round(colMeans(delta.MAP), 3)    # Est
Table.6[i2+6:9, 5] <- round(colMeans(delta.SE), 3)     # SE

## Avg. number of tests:
Table.5[i1+11, 4] <- round( mean(T), 1 )    
Table.5[i1+12, 4] <- round( sd(T), 1 )

#========= Estimation with AT protocol =========#

protocol <- AT
list.zmat <- lapply( rep(N, nsims), array.2dim.data, 
                     design=protocol, Se=Se, Sp=Sp, Yt=Yt )

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
Table.5[i1+1:4, 7] <- round(colMeans(p.Mean), 3)   # Est
Table.5[i1+1:4, 8] <- round(colMeans(p.SE), 4)     # SE

Table.6[i2+1:4, 7] <- round(colMeans(delta.Mean), 3)   # Est
Table.6[i2+1:4, 8] <- round(colMeans(delta.SE), 3)     # SE

## MAP Estimation:
Table.5[i1+6:9, 7] <- round(colMeans(p.MAP), 3)    # Est
Table.5[i1+6:9, 8] <- round(colMeans(p.SE), 4)     # SE

Table.6[i2+6:9, 7] <- round(colMeans(delta.MAP), 3)    # Est
Table.6[i2+6:9, 8] <- round(colMeans(delta.SE), 3)     # SE

## Avg. number of tests:
Table.5[i1+11, 7] <- round( mean(T), 1 )    
Table.5[i1+12, 7] <- round( sd(T), 1 )

}

stopCluster( cl )

## Preparing & saving Table 5
dfr1 <- data.frame( cbind(
     c("","","","","Urine","N=4402","","","","","","", 
	   "",  "","","","","Swab","N=10048","","","","","",""),  
     c("Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T", "", 
	   "Mean", "", "", "", "", "MAP", "", "", "", "", "T.bar", "S.T"),
     c("-", "+", "-", "+", "", "-", "+", "-", "+", "", "", "", "", 
       "-", "+", "-", "+", "", "-", "+", "-", "+", "", "", ""),
     c("-", "-", "+", "+", "", "-", "-", "+", "+", "", "", "", "", 
       "-", "-", "+", "+", "", "-", "-", "+", "+", "", "", "")
	   ))

dfr1 <- data.frame( dfr1, Table.5 )
dfr1 <- sapply(dfr1, as.character)
dfr1[is.na(dfr1)] <- " "
Table_5 <- as.data.frame(dfr1)
## Mean
T <- as.numeric( c( Table_5[11, 5], Table_5[11, 8], Table_5[11, 11] ) )
redc <- round( 100*(T[1] - T[-1])/T[1], 1 )
T.redc <- c( T[1],
             paste( T[2], "(", redc[1], "%", ")", sep="" ),
             paste( T[3], "(", redc[2], "%", ")", sep="" )
            )
Table_5[11, 8] <- T.redc[2]
Table_5[11, 11] <- T.redc[3]

## MAP
T <- as.numeric( c( Table_5[24, 5], Table_5[24, 8], Table_5[24, 11] ) )
redc <- round( 100*(T[1] - T[-1])/T[1], 1 )
T.redc <- c( T[1],
             paste( T[2], "(", redc[1], "%", ")", sep="" ),
             paste( T[3], "(", redc[2], "%", ")", sep="" )
            )
Table_5[24, 8] <- T.redc[2]
Table_5[24, 11] <- T.redc[3]
colnames( Table_5 ) <- c("Stratum", "", "CT", "NG", "Est", "SE", "", "Est", "SE", "", "Est", "SE")
                            				   
saveRDS(Table_5, "./data_analysis/intermediate_results/Table_5.rds")

## Preparing & saving Table 6
dfr2 <- data.frame( cbind(
     c("","","","","Urine","N=4402","","","","","","", 
	   "","","Swab","N=10048","","",""),
     c("", "", "", "Se:(1)1=0.947", "Se:(1)2=0.913", "Sp:(1)1=0.989", "Sp:(1)2=0.993", "", "", "",  
	   "", "", "", "Se:(1)1=0.942", "Se:(1)2=0.992", "Sp:(1)1=0.976", "Sp:(1)2=0.987", "", ""),  
     c("Mean", "", "", "", "", "MAP", "", "", "", "",  
	   "Mean", "", "", "", "", "MAP", "", "", "")))

dfr2 <- data.frame( dfr2, Table.6 )
dfr2 <- sapply(dfr2, as.character)
dfr2[is.na(dfr2)] <- " "

Table_6 <- as.data.frame(dfr2)
colnames( Table_6 ) <- c("Stratum", "Accuracy", "", "Est", "SE", "", "Est", "SE", "", "Est", "SE")
                            
						   
saveRDS(Table_6, "./data_analysis/intermediate_results/Table_6.rds")

