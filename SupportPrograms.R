
#******************************************************#
# Support function to implement the EM algorithm.
# Assay accuracies: known.
#
EmKnownAssayAcr <- function(p0,Z,Yt,N,S,p.pr,emGit,emburn,emmaxit,emtol,tracing){
  dyn.load("./source_code/gbbstwodisgen.dll")  # Call FORTRAN program
  Ytmat <- poolMembTracker(Z,Yt,N,S)
  N.til <- N-4+sum(p.pr)
  Ycol <- ncol(Ytmat)
  SeSp <- Z[ ,4:7]
  Z <- Z[ ,-(4:7)]
  Zrow <- nrow(Z)
  Zcol <- ncol(Z)
  p1 <- p0
  p0 <- p0 + 2*emtol
  s <- 1
  convergence <- 0
  while(max(abs(p1-p0)) > emtol){  # Start iterating
    p0 <- p1
    V <- matrix(0,nrow=N,ncol=4)
    U <- matrix(runif(N*emGit),nrow=N,ncol=emGit)
    p1 <- EmultGibbs(V,p0,Ytmat,Z,N,SeSp,Ycol,Zrow,Zcol,U,emGit,emburn)$res # E step
	p1 <- p1/N.til + (p.pr-1)/N.til  # Updating by prior
	if(s >= emmaxit){
      convergence <- 1
	  break
    }
    s <- s + 1
	if(tracing){
      print(c(s,p1))
	}
  }
  dyn.unload("./source_code/gbbstwodisgen.dll")
  list("prevalence"=p1,"accuracy"=NULL,"convergence"=convergence)
}


#******************************************************#
# Support function to implement the EM algorithm.
# Assay accuracies: unknown.
#
EmUnknownAssayAcr <- function(p0,delta0,Z,Yt,N,S,p.pr,se1.pr,se2.pr,
                     sp1.pr,sp2.pr,emGit,emburn,emmaxit,emtol,tracing){
  dyn.load("./source_code/mapacrtwodgen.dll")  # Call FORTRAN program
  Ytmat <- poolMembTracker(Z,Yt,N,S)
  N.til <- N-4+sum(p.pr)
  Ycol <- ncol(Ytmat)
  SeSp <- Z[ ,4:7]
  Z <- Z[ ,-(4:7)]
  Zrow <- nrow(Z)
  Zcol <- ncol(Z)
  J <- Zrow
  param0 <- c(p0,delta0)
  param1 <- param0
  param0 <- param0 + 2*emtol
  s <- 1
  convergence <- 0
  while(max(abs(param1-param0)) > emtol){  # Start iterating
    param0 <- param1
	p0 <- param0[1:4]
    SeSp[ ,1] <- param0[5]
    SeSp[ ,2] <- param0[6]
    SeSp[ ,3] <- param0[7]
    SeSp[ ,4] <- param0[8]
    V <- matrix(0,nrow=N,ncol=4)
    U <- matrix(runif(N*emGit),nrow=N,ncol=emGit)
	res <- EmultGibbs2(V,p0,Ytmat,Z,N,SeSp,Ycol,Zrow,Zcol,U,emGit,emburn)  # E step
	p1 <- (p.pr-1)/N.til + res$emul/N.til   # Updating by prior
	eztil <- res$eztil
	J.se1 <- sum(se1.pr)-2 + sum(eztil[ ,1])
	J.se2 <- sum(se2.pr)-2 + sum(eztil[ ,2])
	J.sp1 <- sum(sp1.pr)-2 + J-sum(eztil[ ,1])
	J.sp2 <- sum(sp2.pr)-2 + J-sum(eztil[ ,2])
	se1 <- (se1.pr[1]-1 + sum(Z[ ,1]*eztil[ ,1]) )/J.se1
	se2 <- (se2.pr[1]-1 + sum(Z[ ,2]*eztil[ ,2]) )/J.se2
	sp1 <- (sp1.pr[1]-1 + sum( (1-Z[ ,1])*(1-eztil[ ,1]) ) )/J.sp1
	sp2 <- (sp2.pr[1]-1 + sum( (1-Z[ ,2])*(1-eztil[ ,2]) ) )/J.sp2
    param1 <- c(p1,se1,se2,sp1,sp2)
	if(s >= emmaxit){
      convergence <- 1
	  break
    }
    s <- s + 1
	if(tracing){
      print(c(s,param1))
	}
  }
  dyn.unload("./source_code/mapacrtwodgen.dll")
  list("prevalence"=param1[1:4],"accuracy"=param1[5:8],"convergence"=convergence)
}

#******************************************************#
# Support function to approximate expectations using 
# a Gibbs sampler in the E step of the EM algorithm. 
#
EmultGibbs <- function(V,p,Ytmat,Z,N,SeSp,Ycol,Zrow,Zcol,U,emGit,emburn){
  res <- .C("gbbstwodisgen",as.integer(V),p,as.integer(Ytmat),
            as.integer(Z),as.integer(N),SeSp,as.integer(Ycol),
            as.integer(Zrow),as.integer(Zcol),U,as.integer(emGit),
            as.integer(emburn))
  res2 <- colSums( matrix(res[[1]],N,4)/(emGit-emburn) )
  ytil <- res[[3]]
  list("res"=res2,"ytil"=ytil)
}

#******************************************************#
# Support function to approximate expectations using 
# a Gibbs sampler in the E step of the EM algorithm. 
#
EmultGibbs2 <- function(V,p,Ytmat,Z,N,SeSp,Ycol,Zrow,Zcol,U,emGit,emburn){
  Ztil <- matrix(0,Zrow,2)
  res <- .C("mapacrtwodgen",as.integer(V),p,as.integer(Ytmat),
            as.integer(Z),as.integer(Ztil),as.integer(N),SeSp,
			as.integer(Ycol),as.integer(Zrow),as.integer(Zcol),U,
			as.integer(emGit),as.integer(emburn))	
  emul <- colSums( matrix(res[[1]],N,4)/(emGit-emburn) )
  eztil <- matrix(res[[5]],Zrow,2)/(emGit-emburn)
  list("emul"=emul,"eztil"=eztil)
}

#******************************************************#
# Support function to implement the Bayesian methods. 
# Assay accuracies: known.
#
PostKnownAssayAcr <- function(p0,Z,Yt,N,S,p.pr,postGit,tracing){
  library("MCMCpack") # Used to sample from Dirichlet
  dyn.load("./source_code/ytiltwodbayes.dll")  # Call FORTRAN program
  Ytmat <- poolMembTracker(Z,Yt,N,S)
  Ycol <- ncol(Ytmat)
  SeSp <- Z[ ,4:7]
  Z <- Z[ ,-(4:7)]
  Zrow <- nrow(Z)
  Zcol <- ncol(Z)
  p.save <- sesp.save <- matrix(-9,postGit,4)
  p <- p0
  for(g in 1:postGit){  # Start iterating
    U <- runif(N)
    Zt <- matrix(0,Zrow,2)
    res <- .C("ytiltwodbayes",as.integer(Zt),p,as.integer(Ytmat),
            as.integer(Z),as.integer(N),SeSp,as.integer(Ycol),
            as.integer(Zrow),as.integer(Zcol),U)  # Sample Yt from multinomial
    Ytmat <- matrix(res[[3]],N,Ycol)
    a1 <- p.pr[1] + sum((1-Ytmat[ ,1])*(1-Ytmat[ ,2]))
    a2 <- p.pr[2] + sum(Ytmat[ ,1]*(1-Ytmat[ ,2]))
    a3 <- p.pr[3] + sum((1-Ytmat[ ,1])*Ytmat[ ,2])
    a4 <- p.pr[4] + sum(Ytmat[ ,1]*Ytmat[ ,2])
    p <- MCMCpack::rdirichlet(1,c(a1,a2,a3,a4))  # Sample p from Dirichlet
    p.save[g, ] <- p
	if(tracing){
      print(g)
	}
  }
  dyn.unload("./source_code/ytiltwodbayes.dll")
  colnames(p.save) <- c("p00","p10","p01","p11")
  list("prevalence"=p.save,"accuracy"=NULL,"convergence"=0)
}

#******************************************************#
# Support function to implement the Bayesian methods. 
# Assay accuracies: unknown.
#
PostUnknownAssayAcr <- function(p0,delta0,Z,Yt,N,S,p.pr,se1.pr,
                                se2.pr,sp1.pr,sp2.pr,postGit,tracing){
  library("MCMCpack") # Used to sample from Dirichlet
  dyn.load("./source_code/ytiltwodbayes.dll")  # Call FORTRAN program
  Ytmat <- poolMembTracker(Z,Yt,N,S)
  Ycol <- ncol(Ytmat)
  SeSp <- Z[ ,4:7]
  Z <- Z[ ,-(4:7)]
  Zrow <- nrow(Z)
  Zcol <- ncol(Z)
  p.save <- sesp.save <- matrix(-9,postGit,4)
  p <- p0
  for(g in 1:postGit){
    U <- runif(N)
    Zt <- matrix(0,Zrow,2)
    res <- .C("ytiltwodbayes",as.integer(Zt),p,as.integer(Ytmat),
            as.integer(Z),as.integer(N),SeSp,as.integer(Ycol),
            as.integer(Zrow),as.integer(Zcol),U)  # Sample Yt from multinomial
    Ytmat <- matrix(res[[3]],N,Ycol)
    Ztil <-  matrix(res[[1]],Zrow,2)
    a1 <- p.pr[1] + sum((1-Ytmat[ ,1])*(1-Ytmat[ ,2]))
    a2 <- p.pr[2] + sum(Ytmat[ ,1]*(1-Ytmat[ ,2]))
    a3 <- p.pr[3] + sum((1-Ytmat[ ,1])*Ytmat[ ,2])
    a4 <- p.pr[4] + sum(Ytmat[ ,1]*Ytmat[ ,2])
    p <- MCMCpack::rdirichlet(1,c(a1,a2,a3,a4)) # Sample from Dirichlet
    p.save[g, ] <- p
    a.se1 <- se1.pr[1] + sum(Z[ ,1]*Ztil[ ,1])
    b.se1 <- se1.pr[2] + sum((1-Z[ ,1])*Ztil[ ,1])
    se1 <- rbeta(1,a.se1,b.se1)  # Sample from beta
    a.se2 <- se2.pr[1] + sum(Z[ ,2]*Ztil[ ,2])
    b.se2 <- se2.pr[2] + sum((1-Z[ ,2])*Ztil[ ,2])
    se2 <- rbeta(1,a.se2,b.se2)
    a.sp1 <- sp1.pr[1] + sum((1-Z[ ,1])*(1-Ztil[ ,1]))
    b.sp1 <- sp1.pr[2] + sum(Z[ ,1]*(1-Ztil[ ,1]))
    sp1 <- rbeta(1,a.sp1,b.sp1)
    a.sp2 <- sp2.pr[1] + sum((1-Z[ ,2])*(1-Ztil[ ,2]))
    b.sp2 <- sp2.pr[2] + sum(Z[ ,2]*(1-Ztil[ ,2]))
    sp2 <- rbeta(1,a.sp2,b.sp2)
    Se <- c(se1,se2)
    Sp <- c(sp1,sp2)
    SeSp[ ,1] <- se1
    SeSp[ ,2] <- se2
    SeSp[ ,3] <- sp1
    SeSp[ ,4] <- sp2
    sesp.save[g, ] <- c(Se,Sp)
    if(tracing){
      print(g)
    }
  }
  dyn.unload("./source_code/ytiltwodbayes.dll")
  colnames(p.save) <- c("p00","p10","p01","p11")
  colnames(sesp.save) <- c("Se1","Se2","Sp1","Sp2")
  list("prevalence"=p.save,"accuracy"=sesp.save,"convergence"=0)
}

#******************************************************#
# Support function to keep track of the pool members.
#
poolMembTracker <- function(Z,Yt,N,S){
  K <- 2
  ytm <- matrix(-9,N,S)
  tmp <- Z[ ,-(1:(K+1+2*K))]
  vec <- 1:nrow(tmp)
  nz <- ncol(tmp)
  for(d in 1:N){
    tid <- tmp==d
    store <- NULL
    for(i in 1:nz){
      store <- c(store,vec[tid[ ,i]])
    }
    ytm[d,1:length(store)] <- sort(store)
  }
  cbind(Yt,rowSums(ytm>0),ytm)
}


#########################################################
#########################################################


#******************************************************#
# Support function to implement the EM algorithm.
# Assay accuracies: unknown.
# 
EmUnknownAssayAcr_L2 <- function(p0,delta0,Z,Yt,N,S,p.pr,se.pr,sp.pr,
                                 emGit,emburn,emmaxit,emtol){
  dyn.load("./source_code/mapacrtwodgen.dll")  # Call FORTRAN program

  Ytmat <- poolMembTracker2_L2(Z[ ,-(1:8)],Yt,N,S)
  N.til <- N-4+sum(p.pr)
  Ycol <- ncol(Ytmat)
# SeSp <- Z[ ,4:7]
  Lid <- Z[ ,8]
  L <- length(unique(Lid))
  Z <- Z[ ,-(4:8)]
  Zrow <- nrow(Z)
  Zcol <- ncol(Z)
  SeSp <- matrix(-9,Zrow,4)

  param0 <- c(p0,delta0)
  param1 <- param0
  param0 <- param0 + 2*emtol
  s <- 1
  convergence <- 0
  while(max(abs(param1-param0)) > emtol){  # Start iterating
    param0 <- param1
	p0 <- param0[1:4]

	Mid1 <- 1==Lid
    SeSp[Mid1,1] <- param0[5]
    SeSp[Mid1,2] <- param0[6]
    SeSp[Mid1,3] <- param0[7]
    SeSp[Mid1,4] <- param0[8]
	Mid2 <- 2==Lid
    SeSp[Mid2,1] <- param0[9]
    SeSp[Mid2,2] <- param0[10]
    SeSp[Mid2,3] <- param0[11]
    SeSp[Mid2,4] <- param0[12]

    V <- matrix(0,nrow=N,ncol=4)
    U <- matrix(runif(N*emGit),nrow=N,ncol=emGit)
	res <- EmultGibbs2_L2(V,p0,Ytmat,Z,N,SeSp,Ycol,Zrow,Zcol,U,emGit,emburn)  # E step
	p1 <- (p.pr-1)/N.til + res$emul/N.til   # Updating by prior
	eztil <- res$eztil

	Miid1 <- 1==Lid
	J.se1 <- sum(se.pr[1,1:2])-2 + sum(eztil[Miid1 ,1])
	J.se2 <- sum(se.pr[1,3:4])-2 + sum(eztil[Miid1 ,2])
	J.sp1 <- sum(sp.pr[1,1:2])-2 + sum(1-eztil[Miid1 ,1])
	J.sp2 <- sum(sp.pr[1,3:4])-2 + sum(1-eztil[Miid1 ,2])
	se1 <- (se.pr[1,1]-1 + sum(Z[Miid1 ,1]*eztil[Miid1 ,1]) )/J.se1
	se2 <- (se.pr[1,3]-1 + sum(Z[Miid1 ,2]*eztil[Miid1 ,2]) )/J.se2
	sp1 <- (sp.pr[1,1]-1 + sum( (1-Z[Miid1 ,1])*(1-eztil[Miid1 ,1]) ) )/J.sp1
	sp2 <- (sp.pr[1,3]-1 + sum( (1-Z[Miid1 ,2])*(1-eztil[Miid1 ,2]) ) )/J.sp2
	sesp1 <- c(se1,se2,sp1,sp2)

	Miid2 <- 2==Lid
	LJ.se1 <- sum(se.pr[2,1:2])-2 + sum(eztil[Miid2 ,1])
	LJ.se2 <- sum(se.pr[2,3:4])-2 + sum(eztil[Miid2 ,2])
	LJ.sp1 <- sum(sp.pr[2,1:2])-2 + sum(1-eztil[Miid2 ,1])
	LJ.sp2 <- sum(sp.pr[2,3:4])-2 + sum(1-eztil[Miid2 ,2])
	Lse1 <- (se.pr[2,1]-1 + sum(Z[Miid2 ,1]*eztil[Miid2 ,1]) )/LJ.se1
	Lse2 <- (se.pr[2,3]-1 + sum(Z[Miid2 ,2]*eztil[Miid2 ,2]) )/LJ.se2
	Lsp1 <- (sp.pr[2,1]-1 + sum( (1-Z[Miid2 ,1])*(1-eztil[Miid2 ,1]) ) )/LJ.sp1
	Lsp2 <- (sp.pr[2,3]-1 + sum( (1-Z[Miid2 ,2])*(1-eztil[Miid2 ,2]) ) )/LJ.sp2
	sesp2 <- c(Lse1,Lse2,Lsp1,Lsp2)

    param1 <- c(p1,sesp1,sesp2)
	if(s >= emmaxit){
      convergence <- 1
	  break
    }
    s <- s + 1
#    print(c(s,param1))
  }
  dyn.unload("./source_code/mapacrtwodgen.dll")
  list("prevalence"=param1[1:4],"accuracy"=param1[-(1:4)],"convergence"=convergence)
}

#******************************************************#
# Support function to approximate expectations using 
# a Gibbs sampler in the E step of the EM algorithm. 
#
EmultGibbs2_L2 <- function(V,p,Ytmat,Z,N,SeSp,Ycol,Zrow,Zcol,U,emGit,emburn){
  Ztil <- matrix(0,Zrow,2)
  res <- .C("mapacrtwodgen",as.integer(V),p,as.integer(Ytmat),
            as.integer(Z),as.integer(Ztil),as.integer(N),SeSp,
			as.integer(Ycol),as.integer(Zrow),as.integer(Zcol),U,
			as.integer(emGit),as.integer(emburn))	
  emul <- colSums( matrix(res[[1]],N,4)/(emGit-emburn) )
  eztil <- matrix(res[[5]],Zrow,2)/(emGit-emburn)
  list("emul"=emul,"eztil"=eztil)
}

#******************************************************#
# Support function to implement the Bayesian methods. 
# Assay accuracies: unknown.
#
PostUnknownAssayAcr_L2 <- function(p0,delta0,Z,Yt,N,S,p.pr,se.pr,sp.pr,postGit){
  library("MCMCpack") # Used to sample from Dirichlet
  dyn.load("./source_code/ytiltwodbayes.dll")  # Call FORTRAN program

  Ytmat <- poolMembTracker2_L2(Z[ ,-(1:8)],Yt,N,S)
  Ycol <- ncol(Ytmat)
#  SeSp <- Z[ ,4:7]
  Lid <- Z[ ,8]
  L <- length(unique(Lid))
  Z <- Z[ ,-(4:8)]
  Zrow <- nrow(Z)
  Zcol <- ncol(Z)

  param0 <- c(p0,delta0)
  SeSp <- matrix(-9,Zrow,4)

  Mid1 <- 1==Lid
  SeSp[Mid1,1] <- param0[5]
  SeSp[Mid1,2] <- param0[6]
  SeSp[Mid1,3] <- param0[7]
  SeSp[Mid1,4] <- param0[8]

  Mid2 <- 2==Lid
  SeSp[Mid2,1] <- param0[9]
  SeSp[Mid2,2] <- param0[10]
  SeSp[Mid2,3] <- param0[11]
  SeSp[Mid2,4] <- param0[12]
  
  
  p.save <- matrix(-9,postGit,4)
  sesp.save <- matrix(-9,postGit,8)
  p <- p0
  for(g in 1:postGit){
    U <- runif(N)
    Zt <- matrix(0,Zrow,2)
    res <- .C("ytiltwodbayes",as.integer(Zt),p,as.integer(Ytmat),
            as.integer(Z),as.integer(N),SeSp,as.integer(Ycol),
            as.integer(Zrow),as.integer(Zcol),U)  # Sample Yt from multinomial
    Ytmat <- matrix(res[[3]],N,Ycol)
    Ztil <-  matrix(res[[1]],Zrow,2)
    a1 <- p.pr[1] + sum((1-Ytmat[ ,1])*(1-Ytmat[ ,2]))
    a2 <- p.pr[2] + sum(Ytmat[ ,1]*(1-Ytmat[ ,2]))
    a3 <- p.pr[3] + sum((1-Ytmat[ ,1])*Ytmat[ ,2])
    a4 <- p.pr[4] + sum(Ytmat[ ,1]*Ytmat[ ,2])
    p <- rdirichlet(1,c(a1,a2,a3,a4)) # Sample from Dirichlet
    p.save[g, ] <- p

	Miid1 <- 1==Lid
    a.se1 <- se.pr[1,1] + sum(Z[Miid1 ,1]*Ztil[Miid1 ,1])
    b.se1 <- se.pr[1,2] + sum((1-Z[Miid1 ,1])*Ztil[Miid1 ,1])
    se1 <- rbeta(1,a.se1,b.se1)  # Sample from beta

    a.se2 <- se.pr[1,3] + sum(Z[Miid1 ,2]*Ztil[Miid1 ,2])
    b.se2 <- se.pr[1,4] + sum((1-Z[Miid1 ,2])*Ztil[Miid1 ,2])
    se2 <- rbeta(1,a.se2,b.se2)

    a.sp1 <- sp.pr[1,1] + sum((1-Z[Miid1 ,1])*(1-Ztil[Miid1 ,1]))
    b.sp1 <- sp.pr[1,2] + sum(Z[Miid1 ,1]*(1-Ztil[Miid1 ,1]))
    sp1 <- rbeta(1,a.sp1,b.sp1)

    a.sp2 <- sp.pr[1,3] + sum((1-Z[Miid1 ,2])*(1-Ztil[Miid1 ,2]))
    b.sp2 <- sp.pr[1,4] + sum(Z[Miid1 ,2]*(1-Ztil[Miid1 ,2]))
    sp2 <- rbeta(1,a.sp2,b.sp2)

    Se <- c(se1,se2)
    Sp <- c(sp1,sp2)
    SeSp[Miid1, 1] <- se1
    SeSp[Miid1, 2] <- se2
    SeSp[Miid1, 3] <- sp1
    SeSp[Miid1, 4] <- sp2


	Miid2 <- 2==Lid
    La.se1 <- se.pr[2,1] + sum(Z[Miid2 ,1]*Ztil[Miid2 ,1])
    Lb.se1 <- se.pr[2,2] + sum((1-Z[Miid2 ,1])*Ztil[Miid2 ,1])
    Lse1 <- rbeta(1,La.se1,Lb.se1)  # Sample from beta

    La.se2 <- se.pr[2,3] + sum(Z[Miid2 ,2]*Ztil[Miid2 ,2])
    Lb.se2 <- se.pr[2,4] + sum((1-Z[Miid2 ,2])*Ztil[Miid2 ,2])
    Lse2 <- rbeta(1,La.se2,Lb.se2)

    La.sp1 <- sp.pr[2,1] + sum((1-Z[Miid2 ,1])*(1-Ztil[Miid2 ,1]))
    Lb.sp1 <- sp.pr[2,2] + sum(Z[Miid2 ,1]*(1-Ztil[Miid2 ,1]))
    Lsp1 <- rbeta(1,La.sp1,Lb.sp1)

    La.sp2 <- sp.pr[2,3] + sum((1-Z[Miid2 ,2])*(1-Ztil[Miid2 ,2]))
    Lb.sp2 <- sp.pr[2,4] + sum(Z[Miid2 ,2]*(1-Ztil[Miid2 ,2]))
    Lsp2 <- rbeta(1,La.sp2,Lb.sp2)

    LSe <- c(Lse1,Lse2)
    LSp <- c(Lsp1,Lsp2)
    SeSp[Miid2, 1] <- Lse1
    SeSp[Miid2, 2] <- Lse2
    SeSp[Miid2, 3] <- Lsp1
    SeSp[Miid2, 4] <- Lsp2

    sesp.save[g, ] <- c(Se,Sp,LSe,LSp)
#    print(g)
  }
  dyn.unload("./source_code/ytiltwodbayes.dll")
  colnames(p.save) <- c("p00","p10","p01","p11")
  colnames(sesp.save) <- c("Se(1)1","Se(1)2","Sp(1)1","Sp(1)2","Se(2)1","Se(2)2","Sp(2)1","Sp(2)2")
  list("prevalence"=p.save,"accuracy"=sesp.save,"convergence"=0)
}

#******************************************************#
# Support function to keep track of the pool members.
#
poolMembTracker2_L2 <- function(Zsub,Yt,N,S){
  K <- 2
  ytm <- matrix(-9,N,S)
  tmp <- Zsub
  vec <- 1:nrow(tmp)
  nz <- ncol(tmp)
  for(d in 1:N){
    tid <- tmp==d
    store <- NULL
    for(i in 1:nz){
      store <- c(store,vec[tid[ ,i]])
    }
    ytm[d,1:length(store)] <- sort(store)
  }
  cbind(Yt,rowSums(ytm>0),ytm)
}

