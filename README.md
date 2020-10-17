# Multiplex-Group-Testing

This repository contains R programs of the article, "Estimating the prevalence of two or more diseases using outcomes from multiplex group testing." An R function "multDiseaseBayes" is provided to fit the proposed Bayesian estimation method and maximum a posteriori (MAP) estimation method, which can accommodate ANY group testing data involving multiplex assays and provide cost-effective estimates for the prevalence of multiple diseases as well as the assay accuracies (sensitivity and specificity). 

MainProgram.R - file that consists of the main R function combining the R subfunctions and FORTRAN subroutines. Documentation with illustrative (simulation) examples is also provided. 

SupportPrograms.txt - file that consists of standalone, executable R functions.

Note: For efficient execution, compute-intensive parts of the program are written in FORTRAN 95 and called from R through three DLL files, gbbstwodisgen.dll, mapacrtwodgen.dll, and ytiltwodbayes.dll, which work in a 64-bit R package. 


Reference

Warasi, M., McMahan, C., Tebbs, J., and Bilder, C. (2020+). Estimating the prevalence of two or more diseases using outcomes from multiplex group testing. Under review.


################## Simulation examples ##################

## Specify the working directory
setwd(dir = "C:\\programs")

## Import source files
source("SupportPrograms.txt")
source("MultistageHierarchicalData.txt")
source("TwoStageArrayData.txt")

## Examles from Section 5 of the article with one assay, L = 1, are as follows.

## Simulation configuration I:
N <- 5000
p <- c(.95,.02,.02,.01)
Se <- c(.95,.95)
Sp <- c(.99,.99)
design <- c(9,3,1)  # Three-stage hierarchical design

set.seed(123)
out <- hier.alg.data(p,N,design,Se,Sp)
Z <- out$Data
T <- out$T
## MAP estimation with unknown accuracies and flat priors
res <- multDiseaseBayes(p0=c(.90,.06,.03,.01),delta0=c(.95,.95,.98,.98),
                  Z=Z,Yt=matrix(0,N,2),N=N,S=length(design),N0=0,a0=0,
                  b0=0,acr.info=matrix(0,2,4),emGit=15000,emburn=5000,
	                emmaxit=200,emtol=1e-04,method="MAP",accuracy="unknown")

## MAP Results (equivalent to MLE with these flat priors):
# > res
# $prevalence
# [1] 0.95242612 0.01885368 0.01903360 0.00968660
# 
# $accuracy
# [1] 0.9394198 0.9526531 0.9952553 0.9924425
# 
# $convergence
# [1] 0

set.seed(123)
out <- hier.alg.data(p,N,design,Se,Sp)
Z <- out$Data
T <- out$T
## Bayesian estimates with unknown accuracies and flat priors
res <- multDiseaseBayes(p0=c(.90,.06,.03,.01),delta0=c(.95,.95,.98,.98),
                  Z=Z,Yt=matrix(0,N,2),N=N,S=length(design),N0=0,a0=0,
                  b0=0,acr.info=matrix(0,2,4),postGit=15000,
                  method="Bayesian",accuracy="unknown")

## Bayesian results:
burn <- 5000   # burn-in period
colMeans( res$prevalence[-(1:burn), ] )
colMeans( res$accuracy[-(1:burn), ] )
apply(res$prevalence[-(1:burn),],2,sd)
apply(res$accuracy[-(1:burn),],2,sd)
#
# > colMeans( res$prevalence[-(1:burn), ] )
#         p00         p10         p01         p11 
# 0.951712405 0.019132936 0.019277557 0.009877102 
#
# > colMeans( res$accuracy[-(1:burn), ] )
#       Se1       Se2       Sp1       Sp2 
# 0.9361918 0.9493255 0.9943147 0.9916562 
#
# > apply(res$prevalence[-(1:burn),],2,sd)
#         p00         p10         p01         p11 
# 0.003297922 0.002194269 0.002149631 0.001407028 
# 
# > apply(res$accuracy[-(1:burn),],2,sd)
#         Se1         Se2         Sp1         Sp2 
# 0.017049577 0.014261079 0.002965814 0.003327777 

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
                  
