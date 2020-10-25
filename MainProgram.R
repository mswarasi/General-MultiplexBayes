#************************************************************************#
# Description
# -----------
# General-purpose estimation with multiplex assays. This function implements 
# the EM algorithm (MAP estimation) and the posterior sampling algorithm in 
# the article, "Estimating the prevalence of two or more diseases using outcomes 
# from multiplex group testing," by Warasi et al. (2020+). The function works 
# with two infections (i.e., K = 2) and one multiplex assay (i.e., L = 1). 
# See the simulation setting in the article (Section 5).
#
#
# Usage
# -----
# multDiseaseBayes(p0=c(.90,.06,.03,.01),delta0=c(.95,.95,.98,.98),Z,
#                  Yt=matrix(0,N,2),N,S,p.pr=rep(1,4),Se1.pr=c(1,1),
#                  Se2.pr=c(1,1),Sp1.pr=c(1,1),Sp2.pr=c(1,1),postGit=6000,
#                  emGit=6000,emburn=1000,emmaxit=200,emtol=1e-03,
#                  method=c("MAP","Mean"),accuracy=c("unknown","known"))
#
#
# Arguments
# ---------
# p0       : The initial value of p=(p00,p10,p01,p11), an estimate from historical data.  
# delta0   : The initial value of delta, an estimate from the assay product 
#            literature. This is used only when the assay accuracies are unknown.
# Z        : A matrix of the observed group testing data, Z. See the details.
# Yt       : A N by 2 matrix of the individual true binary statuses.
# N        : The number of individuals tested (i.e., sample size).
# S        : The maximum number of times an individual may be tested in pools or individually. 
#            For example, S for a hierarchical algorithm is the number of hierarchical stages 
#            and S is 3 for a two-dimensional array algorithm. 
# p.pr     : Dirichlet prior for p.
# Se1.pr   : Beta prior for Se1.
# Se2.pr   : Beta prior for Se2.
# Sp1.pr   : Beta prior for Sp1.
# Sp2.pr   : Beta prior for Sp2.
# postGit  : The number of Gibbs samples to be drawn from the posterior dist.
# emGit    : The number of Gibbs samples to be used in the EM algorithm.
# emburn   : The initial Gibbs samples to be discarded in the EM algorithm.
# emmaxit  : The maximum number of iterations the EM algorithm can run.
# emtol    : The convergence tolerance used in the EM algorithm.
# method   : The estimation method to be used, either "MAP" or "Mean". Defaults to "MAP".
# accuracy : Whether the assay accuracies are known or unknown. Defaults to "unknown".
#
#
# Details
# -------
# Compute-intensive parts of the program are written in FORTRAN and called from R through 
# three DLL files, gbbstwodisgen.dll, mapacrtwodgen.dll, and ytiltwodbayes.dll, which work 
# in a 64-bit R package. To use the R function, download the DLL's, save them in a  
# folder in the computer, and specify the working directory. 
#
#
# The observed data matrix, Z, needs to be in the manner shown below. For illustration, 
# we show part of a three-stage hierarchical pooling data using design 6:2:1.
#
# > head( Z )
#        Z1 Z2 psz  Se1  Se2  Sp1  Sp2 Indv1 Indv2 Indv3 Indv4 Indv5 Indv6
# Pool:1  0  0   6 0.92 0.95 0.96 0.99     1     2     3     4     5     6
# Pool:2  1  1   6 0.92 0.95 0.96 0.99     7     8     9    10    11    12
# Pool:3  1  0   6 0.92 0.95 0.96 0.99    13    14    15    16    17    18
# Pool:4  1  1   6 0.92 0.95 0.96 0.99    19    20    21    22    23    24
# Pool:5  1  0   6 0.92 0.95 0.96 0.99    25    26    27    28    29    30
# Pool:6  0  0   6 0.92 0.95 0.96 0.99    31    32    33    34    35    36
#
# > tail( Z )
#          Z1 Z2 psz  Se1  Se2  Sp1  Sp2 Indv1 Indv2 Indv3 Indv4 Indv5 Indv6
# Pool:432  0  0   1 0.92 0.95 0.96 0.99   859    -9    -9    -9    -9    -9
# Pool:433  0  1   1 0.92 0.95 0.96 0.99   860    -9    -9    -9    -9    -9
# Pool:434  1  0   1 0.92 0.95 0.96 0.99   871    -9    -9    -9    -9    -9
# Pool:435  0  0   1 0.92 0.95 0.96 0.99   872    -9    -9    -9    -9    -9
# Pool:436  0  0   1 0.92 0.95 0.96 0.99   981    -9    -9    -9    -9    -9
# Pool:437  1  1   1 0.92 0.95 0.96 0.99   982    -9    -9    -9    -9    -9
#
# The pool indices are shown on the row name.
# Columns 1-2 : Pool testing responses for infections 1 and 2.
# Column  3   : The pool size used.
# Columns 4-5 : Sensitivities for infections 1 and 2.
# Columns 6-7 : Specificities for infections 1 and 2.
# Columns 8-  : Indices of the individuals assigned to the pool, starting
#               at column 8. When a pool size is smaller than the maximum
#               pool size used, the remaining entries must be filled with
#               -9 (as shown above) or by any NEGATIVE INTEGER.
#
# Note: The function uses the package "MCMCpack" to sample from Dirichlet.
#
#
# Value
# -----
# prevalence  : An estimate of p, a point estimate for MAP but MCMC sample for Mean.
# accuracy    : An estimate of delta, a point estimate for MAP but MCMC sample for Mean.
# convergence : An indicator, either 0 or 1, where 0 indicates successful completion and 1 
#               indicates the iteration reaches emmaxit. For Mean, convergence is always 0.
#
#
multDiseaseBayes <- function(p0=c(.90,.06,.03,.01),delta0=c(.95,.95,.98,.98),
                      Z,Yt=matrix(0,N,2),N,S,p.pr=rep(1,4),Se1.pr=c(1,1),
                      Se2.pr=c(1,1),Sp1.pr=c(1,1),Sp2.pr=c(1,1),postGit=6000,
                      emGit=6000,emburn=1000,emmaxit=200,emtol=1e-03,
                      method=c("MAP","Mean"),accuracy=c("unknown","known")){
  method <- match.arg(method)
  accuracy <- match.arg(accuracy)
  if(method=="MAP"){
    if(accuracy=="known"){
      res <- EmKnownAssayAcr(p0=p0,Z=Z,Yt=Yt,N=N,S=S,p.pr=p.pr,emGit=emGit,
                             emburn=emburn,emmaxit=emmaxit,emtol=emtol)
    }
    if(accuracy=="unknown"){
	  res <- EmUnknownAssayAcr(p0=p0,delta0=delta0,Z=Z,Yt=Yt,N=N,S=S,p.pr=p.pr,
                               se1.pr=Se1.pr,se2.pr=Se2.pr,sp1.pr=Sp1.pr,sp2.pr=Sp2.pr,
                               emGit=emGit,emburn=emburn,emmaxit=emmaxit,emtol=emtol)
	}
  }
  if(method=="Mean"){
    if(accuracy=="known"){
	  res <- PostKnownAssayAcr(p0=p0,Z=Z,Yt=Yt,N=N,S=S,p.pr=p.pr,postGit=postGit)
	}
	if(accuracy=="unknown"){
	  res <- PostUnknownAssayAcr(p0=p0,delta0=delta0,Z=Z,Yt=Yt,N=N,S=S,p.pr=p.pr,
                                 se1.pr=Se1.pr,se2.pr=Se2.pr,sp1.pr=Sp1.pr,
                                 sp2.pr=Sp2.pr,postGit=postGit)
	}
  }
  return(res)
}


##################################################

## Specify the working directory
# setwd(dir = "C:\\programs")
setwd("C:/Users/User/Desktop/GeneralMultilexPooling_102520--GitHub")

## Import source files
source("SupportPrograms.txt")
source("MultistageHierarchicalData.txt")
source("TwoStageArrayData.txt")

##################################################

# Examples from Section 5 of the article
# --------------------------------------

## Simulation configuration I:
N <- 5000
p <- c(.95,.02,.02,.01)
Se <- c(.95,.95)
Sp <- c(.99,.99)
design <- c(9,3,1)  # Three-stage hierarchical design

set.seed(123)
out <- hier.alg.data(p,N,design,Se,Sp)
Z <- out$Data     # the pooling data generated
T <- out$T        # the number of tests expended

## MAP estimation with unknown accuracies and flat priors:
res <- multDiseaseBayes(p0=c(.90,.06,.03,.01),delta0=c(.95,.95,.98,.98),
                  Z=Z,Yt=matrix(0,N,2),N=N,S=length(design),p.pr=rep(1,4),
                  Se1.pr=c(1,1),Se2.pr=c(1,1),Sp1.pr=c(1,1),Sp2.pr=c(1,1),
                  emGit=12000,emburn=2000,emmaxit=200,emtol=1e-03,
                  method="MAP",accuracy="unknown")
	  

## MAP Results (equivalent to MLE with these flat priors):
# > res
# $prevalence
# [1] 0.95234086 0.01895446 0.01899814 0.00970654

# $accuracy
# [1] 0.9379753 0.9526409 0.9955003 0.9923724

# $convergence
# [1] 0


set.seed(123)
out <- hier.alg.data(p,N,design,Se,Sp)
Z <- out$Data
T <- out$T
## Mean estimation with unknown accuracies and flat priors:
res <- multDiseaseBayes(p0=c(.90,.06,.03,.01),delta0=c(.95,.95,.98,.98),
                  Z=Z,Yt=matrix(0,N,2),N=N,S=length(design),p.pr=rep(1,4),
                  Se1.pr=c(1,1),Se2.pr=c(1,1),Sp1.pr=c(1,1),Sp2.pr=c(1,1),
                  postGit=12000,method="Mean",accuracy="unknown")

## Mean estimation results:
burn <- 2000            # burn-in period
pick <- seq(1,10000,5)  # thinning

colMeans( res$prevalence[-(1:burn), ][pick, ] )
colMeans( res$accuracy[-(1:burn), ][pick, ] )

apply( res$prevalence[-(1:burn), ][pick, ], 2, sd  )
apply( res$accuracy[-(1:burn), ][pick, ], 2, sd  )

# > colMeans( res$prevalence[-(1:burn), ][pick, ] )
#         p00         p10         p01         p11 
# 0.951657890 0.019214570 0.019251822 0.009875718 

# > colMeans( res$accuracy[-(1:burn), ][pick, ] )
#       Se1       Se2       Sp1       Sp2 
# 0.9362142 0.9491854 0.9944378 0.9915996 

# > apply( res$prevalence[-(1:burn), ][pick, ], 2, sd  )
#         p00         p10         p01         p11 
# 0.003304783 0.002241015 0.002123925 0.001405953 

# > apply( res$accuracy[-(1:burn), ][pick, ], 2, sd  )
#         Se1         Se2         Sp1         Sp2 
# 0.016963091 0.014256603 0.002898824 0.003337798


## Also try other designs and estimation settings
# design <- c(5,1)        # Two-stage hierarchical
# design <- c(18,6,3,1)   # Four-stage hierarchical

## For a two-dimensional array of dimensions 11x11, try
design <- c(11,11,1)    
out <- array.2dim.data(p,N,design,Se,Sp)
Z <- out$Data
T <- out$T

## MAP with known accuracies and flat Dirichlet
res <- multDiseaseBayes(p0=c(.90,.06,.03,.01),Z=Z,Yt=matrix(0,N,2),N=N,
                  S=length(design),N0=0,a0=0,emGit=15000,emburn=5000,
	          emmaxit=100,emtol=1e-03,method="MAP",accuracy="known")

## Bayesian with known accuracies and flat Dirichlet
res <- multDiseaseBayes(p0=c(.90,.06,.03,.01),Z=Z,Yt=matrix(0,N,2),N=N,
                  S=length(design),N0=0,a0=0,postGit=15000,method="Bayesian",
                  accuracy="known")


