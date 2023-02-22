# General-MultiplexBayes

The R programs in this repository correspond to the article titled "Estimating the prevalence of two or more diseases using outcomes from multiplex group testing," which presents simulation and data analysis results for estimating the prevalence of two infections and the probabilities of assay accuracy (sensitivities and specificities). We offer two documented R functions, namely, mult.gt.bayes() and mult.gt.bayes_L2(), each of which can be utilized to implement the posterior sampling algorithm and the EM algorithm described in Warasi et al. (2023+). The first function can fit the model with L=1 assay, while the second one can be used with L=2 assays.

To reproduce all results of the article, use the zip folder "bimj.202200270-sup-code-and-data.zip," which contains the code and data necessary to replicate the study. Note that the code is designed for use on a 64-bit R installation on a Windows machine. The computationally-intensive portions of the code are implemented in FORTRAN and called into R using the .C function. The code is illustrated below using synthetic data, which are comparable to the chlamydia and gonorrhea testing data obtained from the State Hygienic Laboratory at the University of Iowa.


Reference

Warasi, M., Tebbs, J., McMahan, C., and Bilder, C. (2023+). Estimating the prevalence of two or more diseases using outcomes from multiplex group testing. Under revision at Biometrical Journal.


##################################################

To illustrate the code, we use synthetic data for chlamydia and gonorrhea that were observed using a three-stage hierarchical (H3) protocol. To begin, you should read in the observed data, import the source code, and provide details about the pool sizes.

#### Set the working directory:
setwd(dir = "C:\\programs")

#### Reading in the data:
Z <- readRDS( "Urine_GT_H3.csv" )

#### Load the code and package:
source("MainFunctions.R")

source("SupportPrograms.R")

library(MCMCpack)

#### We simulate data for the urine stratum, where the number of specimens tested is N = 4402: 
N <- 4402

#### Pool sizes for the H3 protocol:
protocol <- c(6, 2, 1)

In order for the R functions to work properly, it is important that the observed test responses and accompanying information (such as pool size and pool member IDs) are structured in a specific manner. This is necessary to extract the required information for posterior sampling and the EM algorithm. The input data must be a matrix object called "Z," with the first two columns containing test responses for diseases 1 and 2, the third column containing pool sizes, columns 4-5 containing assay sensitivities for diseases 1 and 2, columns 6-7 containing specificities for diseases 1 and 2, and column 8 onwards containing identification numbers for the individuals assigned to each pool. Any missing entries must be filled with -9 or any negative number. Below is an example of what a portion of the data should look like:

> head( Z )

       Z1 Z2 psz   Se1   Se2   Sp1   Sp2 Indv1 Indv2 Indv3 Indv4 Indv5 Indv6
       
Pool:1  0  0   6 0.947 0.947 0.989 0.989     1     2     3     4     5     6

Pool:2  1  0   6 0.947 0.947 0.989 0.989     7     8     9    10    11    12

Pool:3  0  0   6 0.947 0.947 0.989 0.989    13    14    15    16    17    18

Pool:4  0  1   6 0.947 0.947 0.989 0.989    19    20    21    22    23    24

Pool:5  1  0   6 0.947 0.947 0.989 0.989    25    26    27    28    29    30

Pool:6  0  0   6 0.947 0.947 0.989 0.989    31    32    33    34    35    36

#### Specify an initial value for p and delta. These values can be specified from historical, pilot study, or any contextual estimate.  
p0 <- c(.92, .05, .02, .01)

delta0 <- c(.96, .96, .98, .98)

#### Specify the prior hyperparameters:
p.pr <- c(1, 1, 1, 1)

Se1.pr <- c(1, 1)

Se2.pr <- c(1, 1)

Sp1.pr <- c(1, 1)

Sp2.pr <- c(1, 1)


#### Specify the values needed for the EM algorithm:

G <- 12000                # number of Gibbs iterates

burn <- 2000              # burn-in period

pick <- seq(1, 10000, 5)  # keeping every 5th

epsilon <- 0.001          # convergence tolerance

emmaxit <- 200            # max number of iterations


##### We illustrate the code with both known and unknown values of the assay accuracy probabilities as follows. 

### Case I: The assay accuracy probabilities are KNOWN

set.seed( 123 )

#### POSTERIOR SAMPLING 
post.samp <- mult.gt.bayes( Z=Z, p0=p0, N=N, S=length(protocol), 
                            p.pr=p.pr, postGit=G, method="Bayes", 
                            accuracy="known")

##### Estimated mean and standard deviation of the posterior distribution for p=(p00, p10, p01, p11):
> round( colMeans(post.samp$prevalence[-(1:burn), ][pick, ]), 4)
   
0.9090 0.0804 0.0061 0.0046 

> round( apply(post.samp$prevalence[-(1:burn), ][pick, ], 2, sd), 4)
   
0.0046 0.0044 0.0012 0.0011 


#### MAP ESTIMATION
res.map <- mult.gt.bayes( Z=Z, p0=p0, N=N, S=length(protocol),
                          p.pr=p.pr, emGit=G, emburn=burn, 
                          method="MAP", accuracy="known" )

##### MAP estimation results for p=(p00, p10, p01, p11):
> round( res.map$prevalence, 4)
        
[1,] 0.9097 0.0801 0.0059 0.0044


##### Number of tests expended:

nrow( Z )

[1] 2373



### Case II: The assay accuracy probabilities are UNKNOWN

In this particular case, the vector of assay accuracy probabilities is assumed to be unknown. In order to accommodate this scenario, we remove the assay sensitivities and specificities from columns 4-7 of the data input matrix Z, and fill in any missing values with -9. The revised data input is shown below.

Z[ ,4:7] <- -9

> head( Z )

       Z1 Z2 psz Se1 Se2 Sp1 Sp2 Indv1 Indv2 Indv3 Indv4 Indv5 Indv6
       
Pool:1  0  0   6  -9  -9  -9  -9     1     2     3     4     5     6

Pool:2  1  0   6  -9  -9  -9  -9     7     8     9    10    11    12

Pool:3  0  0   6  -9  -9  -9  -9    13    14    15    16    17    18

Pool:4  0  1   6  -9  -9  -9  -9    19    20    21    22    23    24

Pool:5  1  0   6  -9  -9  -9  -9    25    26    27    28    29    30

Pool:6  0  0   6  -9  -9  -9  -9    31    32    33    34    35    36



#### POSTERIOR SAMPLING 
post.samp <- mult.gt.bayes( Z=Z, p0=p0, delta0=delta0, N=N, 
                 S=length(protocol), p.pr=p.pr, Se1.pr=Se1.pr, Se2.pr=Se2.pr,
                 Sp1.pr=Sp1.pr, Sp2.pr=Sp2.pr, postGit=G, method="Bayes", 
                 accuracy="unknown" )
                 
                 
##### Estimated mean and standard deviation of the posterior distribution for p=(p00, p10, p01, p11):

> round( colMeans(post.samp$prevalence[-(1:burn), ][pick, ]), 4)

0.9064 0.0831 0.0059 0.0046 

> round( apply(post.samp$prevalence[-(1:burn), ][pick, ], 2, sd), 4)

0.0050 0.0047 0.0012 0.0010 


##### Estimated mean and standard deviation of the posterior distribution for delta=(Se1, Se2, Sp1, Sp2):

> round( colMeans(post.samp$accuracy[-(1:burn), ][pick, ]), 4)

0.9367 0.9686 0.9956 0.9890 

> round( apply(post.samp$accuracy[-(1:burn), ][pick, ], 2, sd), 4)

0.0111 0.0193 0.0035 0.0023 


#### MAP ESTIMATION
res.map <-  mult.gt.bayes( Z=Z, p0=p0, delta0=delta0, N=N, 
                 S=length(protocol), p.pr=p.pr, Se1.pr=Se1.pr, Se2.pr=Se2.pr, 
                 Sp1.pr=Sp1.pr, Sp2.pr=Sp2.pr, emGit=G, emburn=burn, 
                 method="MAP", accuracy="unknown" )


##### MAP estimation results for p=(p00, p10, p01, p11):

> round( res.map$prevalence, 4)

0.9075 0.0825 0.0056 0.0043 

##### MAP estimation results for delta=(Se1, Se2, Sp1, Sp2):

> round( res.map$accuracy, 4)

0.9392 0.9802 0.9961 0.9894 


#### Note: The example shown above is provided in the R script Example.R

