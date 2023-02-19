
#************************************************************************#
# Description
# -----------
# This function can implement the proposed EM algorithm and the posterior 
# sampling algorithm with L=1 assay. It can be used in Windows with 64-bit R.
#
#
# Usage
# -----
# mult.gt.bayes(Z, p0=c(.90,.06,.03,.01), delta0=c(.95,.95,.98,.98),
#               N, S, p.pr=rep(1,4), Se1.pr=c(1,1),
#               Se2.pr=c(1,1), Sp1.pr=c(1,1), Sp2.pr=c(1,1), postGit=6000,
#               emGit=6000, emburn=1000, emmaxit=200, emtol=1e-03,
#               method=c("MAP", "Bayes"), accuracy=c("unknown", "known"))
#
#
# Arguments
# ---------
# Z        : A matrix of the observed group testing data, Z. See the details.
# p0       : Initial value of p=(p00, p10, p01, p11), an estimate from historical data.  
# delta0   : Initial value (4 x 1 vector) of delta, an estimate from the assay product 
#            literature. This is used only when the assay accuracies are unknown.
# N        : Number of individuals tested (i.e., sample size).
# S        : Maximum number of times an individual may be tested in pools or individually. 
#            For example, S for a hierarchical algorithm is the number of hierarchical stages 
#            and S is 3 for a two-dimensional array algorithm. 
# p.pr     : A vector of 4x1 Dirichlet hyperparameters for p.
# Se1.pr   : A vector of 2x1 beta hyperparameters for Se1.
# Se2.pr   : A vector of 2x1 beta hyperparameters for Se2.
# Sp1.pr   : A vector of 2x1 beta hyperparameters for Sp1.
# Sp2.pr   : A vector of 2x1 beta hyperparameters for Sp2.
# postGit  : Number of Gibbs samples to be drawn from the posterior dist.
# emGit    : Number of Gibbs samples to be used in the EM algorithm.
# emburn   : Nnitial Gibbs samples to be discarded in the EM algorithm.
# emmaxit  : Maximum number of iterations the EM algorithm can run.
# emtol    : Convergence tolerance used in the EM algorithm.
# method   : Estimation method to be used, either "MAP" or "Bayes". Defaults to "MAP".
# accuracy : Whether the assay accuracies are known or unknown. Defaults to "unknown".
#
#
# Details
# -------
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
# prevalence  : An estimate of p, a point estimate for MAP but MCMC sample for Bayes.
# accuracy    : An estimate of delta, a point estimate for MAP but MCMC sample for Bayes.
# convergence : An indicator, either 0 or 1, where 0 indicates successful completion and 1 
#               indicates the iteration reaches emmaxit. For Bayes, convergence is always 0.
#
#
mult.gt.bayes <- function(Z, p0=c(.92,.05,.02,.01), delta0=c(.96, .96, .98, .98),
                      N, S, p.pr=rep(1, 4), Se1.pr=c(1, 1),
                      Se2.pr=c(1, 1), Sp1.pr=c(1, 1), Sp2.pr=c(1, 1), postGit=6000,
                      emGit=6000, emburn=1000, emmaxit=200, emtol=1e-03,
                      method=c("MAP", "Bayes"), accuracy=c("unknown", "known")){

  Ymul <- rmultinom(N,1,p0)
  Yt <- cbind(Ymul[2, ] + Ymul[4, ], Ymul[3, ] + Ymul[4, ])

  method <- match.arg(method)
  accuracy <- match.arg(accuracy)
  if(method=="MAP"){
    if(accuracy=="known"){
      res <- EmKnownAssayAcr(p0=p0,Z=Z,Yt=Yt,N=N,S=S,p.pr=p.pr,emGit=emGit,
                             emburn=emburn,emmaxit=emmaxit,emtol=emtol,tracing=FALSE)
    }
    if(accuracy=="unknown"){
	  res <- EmUnknownAssayAcr(p0=p0,delta0=delta0,Z=Z,Yt=Yt,N=N,S=S,p.pr=p.pr,
                               se1.pr=Se1.pr,se2.pr=Se2.pr,sp1.pr=Sp1.pr,sp2.pr=Sp2.pr,
                               emGit=emGit,emburn=emburn,emmaxit=emmaxit,emtol=emtol,tracing=FALSE)
	}
  }
  if(method=="Bayes"){
    if(accuracy=="known"){
	  res <- PostKnownAssayAcr(p0=p0,Z=Z,Yt=Yt,N=N,S=S,p.pr=p.pr,postGit=postGit,tracing=FALSE)
	}
	if(accuracy=="unknown"){
	  res <- PostUnknownAssayAcr(p0=p0,delta0=delta0,Z=Z,Yt=Yt,N=N,S=S,p.pr=p.pr,
                                 se1.pr=Se1.pr,se2.pr=Se2.pr,sp1.pr=Sp1.pr,
                                 sp2.pr=Sp2.pr,postGit=postGit,tracing=FALSE)
	}
  }
  return(res)
}


#************************************************************************#
# Description
# -----------
# This function can implement the proposed EM algorithm and the posterior 
# sampling algorithm with L=2 assays. It can be used in Windows with 64-bit R.
# This function can be used when delta is assumed to be unknown. When delta 
# is known, function mult.gt.bayes() shown above can be used.
#
#
# Usage
# -----
# mult.gt.bayes_L2(Z, p0=c(.92,.05,.02,.01), delta0=c(rep(.96,2),rep(.98,2),rep(.96,2),rep(.98,2)),
#                  N, S, p.pr=rep(1,4), se.pr=matrix(1,2,4), sp.pr=matrix(1,2,4),
#                  postGit=6000, emGit=6000, emburn=1000, emmaxit=100, emtol=1e-03,
#                  method=c("MAP","Bayes"))
#							 
#
# Arguments
# ---------
# Z        : A matrix of the observed group testing data, Z. See the details above.
# p0       : Initial value of p=(p00, p10, p01, p11), an estimate from historical data.  
# delta0   : Initial value (8 x 1 vector) of delta, an estimate from the assay product literature.
# N        : Number of individuals tested (i.e., sample size).
# S        : Maximum number of times an individual may be tested in pools or individually. 
#            For example, S for a hierarchical algorithm is the number of hierarchical stages 
#            and S is 3 for a two-dimensional array algorithm. 
# p.pr     : A vector of 4 x 1 Dirichlet hyperparameters for p.
#
# se.pr    : A 2 x 4 matrix of beta hyperparameters for Se1 and Se2. Columns 1-2 contain
#            the hyperparameters of Se1 and Columns 3-4 the hyperparameters of Se2.
#            Row 1 is for the first assay and row 2 is for the second assay.
#
# sp.pr    : A 2 x 4 matrix of beta hyperparameters for Sp1 and Sp2. Columns 1-2 contain
#            the hyperparameters of Sp1 and Columns 3-4 the hyperparameters of Sp2.
#            Row 1 is for the first assay and row 2 is for the second assay.
#
# postGit  : Number of Gibbs samples to be drawn from the posterior dist.
# emGit    : Number of Gibbs samples to be used in the EM algorithm.
# emburn   : Nnitial Gibbs samples to be discarded in the EM algorithm.
# emmaxit  : Maximum number of iterations the EM algorithm can run.
# emtol    : Convergence tolerance used in the EM algorithm.
# method   : Estimation method to be used, either "MAP" or "Bayes". Defaults to "MAP".
#
#
mult.gt.bayes_L2 <- function(Z, p0=c(.92,.05,.02,.01), 
                             delta0=c(rep(.96,2),rep(.98,2),rep(.96,2),rep(.98,2)),
                             N, S, p.pr=rep(1,4), se.pr=matrix(1,2,4), sp.pr=matrix(1,2,4),
                             postGit=6000, emGit=6000, emburn=1000, emmaxit=100, emtol=1e-03,
                             method=c("MAP", "Bayes")){

  Ymul <- rmultinom(N,1,p0)
  Yt <- cbind(Ymul[2, ] + Ymul[4, ], Ymul[3, ] + Ymul[4, ])

  method <- match.arg(method)
  if(method=="MAP"){
    res <- EmUnknownAssayAcr_L2(p0,delta0,Z,Yt,N,S,p.pr,se.pr,sp.pr,emGit,emburn,emmaxit,emtol)
  }
  if(method=="Bayes"){
    res <- PostUnknownAssayAcr_L2(p0,delta0,Z,Yt,N,S,p.pr,se.pr,sp.pr,postGit)
  }
  return(res)
}

