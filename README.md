# Gen-MultiplexBayes

This repository contains R programs of the article, "Estimating the prevalence of two or more diseases using outcomes from multiplex group testing." An R function "multDiseaseBayes" is provided to fit the proposed Bayesian estimation method and maximum a posteriori (MAP) estimation method, which can accommodate ANY group testing data involving multiplex assays and provide cost-effective estimates for the prevalence of multiple diseases as well as the assay accuracies (sensitivity and specificity). 

Files uploaded:

1. MainProgram.R - R file consisting the main R function that implements the proposed estimation techniques. This combines the subfunctions written in R and FORTRAN. The documentation is also provided.

2. SupportPrograms.txt - This file consists of the supporting R programs.

3. MultistageHierarchicalData.txt - This file provides an R function that simulates multiplex hierarchical group testing data with ANY number of hierarchical stages. This was used for the data simulation in Section 5 of the article. 

4. TwoStageArrayData.txt - This file has an R function that simulates multiplex two-stage array testing data. This was used for the data simulation in Section 5 of the article

5. gbbstwodisgen.dll, mapacrtwodgen.dll, ytiltwodbayes.dll - These DLL files have compiled FORTRAN subroutines for efficient execution of the Gibbs samplers used in the POSTERIOR SAMPLING ALGORITHM as well as in the MAP EM ALGORITHM. The DLL files will work with a 64-bit R package. 


Reference

Warasi, M., McMahan, C., Tebbs, J., and Bilder, C. (2021+). Estimating the prevalence of two or more diseases using outcomes from multiplex group testing. Under review.



###################### SIMULATION EXAMPLES ######################

## Usage

multDiseaseBayes(p0=c(.90,.06,.03,.01),delta0=c(.95,.95,.98,.98),
                      Z,N,S,p.pr=rep(1,4),Se1.pr=c(1,1),
                      Se2.pr=c(1,1),Sp1.pr=c(1,1),Sp2.pr=c(1,1),postGit=6000,
                      emGit=6000,emburn=1000,emmaxit=200,emtol=1e-03,
                      method=c("MAP","Mean"),accuracy=c("unknown","known"))


## Download and save the files in a folder, and specify the directory:
setwd(dir = "C:\\programs")

## Import source files
source("SupportPrograms.txt")

source("MultistageHierarchicalData.txt")

source("TwoStageArrayData.txt")


## Simulation configuration I with one assay (L=1):
N <- 5000

p <- c(.95,.02,.02,.01)

Se <- c(.95,.95)

Sp <- c(.99,.99)

## Example 1: Hierarchical Protocol. MAP with unknown accuracies and flat priors:

set.seed(123)

design <- c(9,3,1)  # Three-stage hierarchical design

out <- hier.alg.data(p,N,design,Se,Sp)

Z <- out$Data   # Simulated pooling data

T <- out$T      # The number of tests expended in the simulation

res <- multDiseaseBayes(p0=c(.90,.06,.03,.01),delta0=c(.95,.95,.98,.98),
                  Z=Z,Yt=matrix(0,N,2),N=N,S=length(design),p.pr=rep(1,4),
                  Se1.pr=c(1,1),Se2.pr=c(1,1),Sp1.pr=c(1,1),Sp2.pr=c(1,1),
                  emGit=12000,emburn=2000,emmaxit=200,emtol=1e-03,
                  method="MAP",accuracy="unknown")

> res

# $prevalence

[1] 0.95235454 0.01894246 0.01900132 0.00970168

# $accuracy

[1] 0.9382384 0.9527329 0.9954653 0.9923734

# $convergence

[1] 0


## Example 2: Hierarchical Protocol. MAP with unknown accuracies and flat priors:

design <- c(9,3,1)  # Three-stage hierarchical design

out <- hier.alg.data(p,N,design,Se,Sp)

Z <- out$Data   # Simulated pooling data

res <- multDiseaseBayes(p0=c(.90,.06,.03,.01),delta0=c(.95,.95,.98,.98),
                  Z=Z,Yt=matrix(0,N,2),N=N,S=length(design),p.pr=rep(1,4),
                  Se1.pr=c(1,1),Se2.pr=c(1,1),Sp1.pr=c(1,1),Sp2.pr=c(1,1),
                  postGit=12000,method="Mean",accuracy="unknown")

burn <- 2000            # burn-in period

pick <- seq(1,10000,5)  # thinning

colMeans( res$prevalence[-(1:burn), ][pick, ] )

colMeans( res$accuracy[-(1:burn), ][pick, ] )

apply( res$prevalence[-(1:burn), ][pick, ], 2, sd  )

apply( res$accuracy[-(1:burn), ][pick, ], 2, sd  )


> colMeans( res$prevalence[-(1:burn), ][pick, ] )

         p00         p10         p01         p11 
	 
 0.951657890 0.019214570 0.019251822 0.009875718 


> colMeans( res$accuracy[-(1:burn), ][pick, ] )

     Se1       Se2       Sp1       Sp2 
      
0.9362142 0.9491854 0.9944378 0.9915996 


> apply( res$prevalence[-(1:burn), ][pick, ], 2, sd  )

     p00         p10         p01         p11 

0.003304783 0.002241015 0.002123925 0.001405953 


> apply( res$accuracy[-(1:burn), ][pick, ], 2, sd  )

     Se1         Se2         Sp1         Sp2 

0.016963091 0.014256603 0.002898824 0.003337798


## Try other hierarchical protocols:

design <- c(5,1)        # Two-stage hierarchical

design <- c(18,6,3,1)   # Four-stage hierarchical


## Example 3: Array Protocol. MAP with known accuracies and flat Dirichlet:

design <- c(11,11,1)    

out <- array.2dim.data(p,N,design,Se,Sp)

Z <- out$Data

res <- multDiseaseBayes(p0=c(.90,.06,.03,.01),Z=Z,Yt=matrix(0,N,2),
                        N=N,S=length(design),p.pr=rep(1,4),
                        emGit=12000,emburn=2000,emmaxit=200,emtol=1e-03,
                        method="MAP",accuracy="known")

> res

# $prevalence
[1] 0.95495022 0.01767218 0.01816344 0.00921416

# $accuracy
NULL

# $convergence
[1] 0
						

## Example 4: Array Protocol. Mean with known accuracies and flat Dirichlet:

design <- c(11,11,1)    

out <- array.2dim.data(p,N,design,Se,Sp)

Z <- out$Data

res <- multDiseaseBayes(p0=c(.90,.06,.03,.01),delta0=c(.95,.95,.98,.98),
                  Z=Z,Yt=matrix(0,N,2),N=N,S=length(design),p.pr=rep(1,4),
                  postGit=12000,method="Mean",accuracy="known")

burn <- 2000            # burn-in period

pick <- seq(1,10000,5)  # thinning

> colMeans( res$prevalence[-(1:burn), ][pick, ] )

     p00         p10         p01         p11 
     
0.954443997 0.017813911 0.018389059 0.009353033 


> apply( res$prevalence[-(1:burn), ][pick, ], 2, sd  )

     p00         p10         p01         p11 
     
0.002954447 0.001928380 0.001923174 0.001411674 


