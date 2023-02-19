
#******************************************************#
# Main function for simulating hierarchical pooling data with L=1 assay
#
hier.alg.data <- function(N,p=NULL,design,Se,Sp,Yt=NULL){
  J <- 2
  c.s <- design
  S <- length(c.s)
  if(S == 1){
    return( mpt.data(p,N,design,Se,Sp,Yt) )
  }
  if(S > 1){
    return( multistage.data(p,N,design,Se,Sp,Yt) )
  }
}

#******************************************************#
# Support function for hier.alg.data()
#
multistage.data <- function(p,N,design,Se,Sp,Yt){
  J <- 2
  c.s <- design
  S <- length(c.s)

  id <- sample( 1:N, N, replace=FALSE )
  Yt <- Yt[ id, ]

  if(!is.null(Yt)){
    Ytil1 <- Yt[ ,1]
    Ytil2 <- Yt[ ,2]
  }
  else{
    Ymul <- rmultinom(N,1,p)
    Ytil1 <- Ymul[2, ] + Ymul[4, ]
    Ytil2 <- Ymul[3, ] + Ymul[4, ]
  } 

  M <- floor(N/c.s[1])
  Rem <- N-M*c.s[1]
  psz <- n.div <- list()
  psz[[1]] <- rep(c.s[1],M)
  n.div[[1]] <- rep(1,length(psz[[1]]))
  n.sm <- matrix(-9,length(psz[[1]]),S)
  n.sm[ ,1] <- 1
  if(S > 1){
    for( s in 1:(S-1) ){
	  store <- tmp <- NULL
	  for(i in 1:length(psz[[s]])){
	    temp <- rep( c.s[s+1], floor(psz[[s]][i]/c.s[s+1]) )
  	    store <- c(store,temp)
	    tmp <- c(tmp,length(temp))
	  }
	  psz[[s+1]] <- store
	  n.div[[s+1]] <- tmp
    }
    vec <- rep(1,length(psz[[1]]))
    for(s in 1:(S-1) ){
	  id0 <- cumsum(c(0,vec))
	  for(i in 1:length(psz[[1]])){
        vec[i] <- sum(n.div[[s+1]][(id0[i]+1):id0[i+1]])
      }
      n.sm[ ,s+1] <- vec
    }
  }
  T <- 0
  Zmat <- NULL
  cc <- cumsum(c(0,colSums(n.sm)))
  id <- cumsum(c(0,psz[[1]]))
  pl.res <- matrix(-9,cc[S+1],2)
  for(m in 1:M){
    mid <- (id[m]+1):id[m+1]
    prob1 <- ifelse(sum(Ytil1[mid])>0,Se[1],1-Sp[1])
    prob2 <- ifelse(sum(Ytil2[mid])>0,Se[2],1-Sp[2])
    z1 <- rbinom(1,1,prob1)
    z2 <- rbinom(1,1,prob2)
    pl.res[m,1] <- z1
    pl.res[m,2] <- z2
    Zmat <- rbind(Zmat,c(z1,z2,length(mid),Se[1],Se[2],Sp[1],Sp[2],mid))
    T <- T + 1
  }
  if( S > 1){
    for(s in 2:S){
      Z1 <- pl.res[(cc[s-1]+1):cc[s],1]
      Z2 <- pl.res[(cc[s-1]+1):cc[s],2]
      cid <- cumsum(c(0,psz[[s]]))
      cn <- cumsum(c(0,n.div[[s]]))
      tmp1 <- tmp2 <- NULL
      for(d in 1:length(psz[[s-1]])){
        tmp3 <- tmp4 <- NULL
        if(Z1[d]==0 & Z2[d]==0){
          tmp3 <- tmp4 <- rep(0,length((cn[d]+1):cn[d+1]))
        }
        if(psz[[s-1]][d]==1){
          tmp3 <- Z1[d]
          tmp4 <- Z2[d]
        }
        if(psz[[s-1]][d]>1){
          if(Z1[d]+Z2[d] > 0){
            for(i in (cn[d]+1):cn[d+1]){
              crng <- (cid[i]+1):cid[i+1]
              prob1 <- ifelse(sum(Ytil1[crng])>0,Se[1],1-Sp[1])
              prob2 <- ifelse(sum(Ytil2[crng])>0,Se[2],1-Sp[2])
              ztp1 <- rbinom(1,1,prob1)
              ztp2 <- rbinom(1,1,prob2)
              tmp3 <- c(tmp3,ztp1)
              tmp4 <- c(tmp4,ztp2)
              fill1 <- rep(-9,psz[[1]][1]-length(crng))
              Zmat <- rbind(Zmat,c(ztp1,ztp2,length(crng),Se[1],
                            Se[2],Sp[1],Sp[2],crng,fill1))
              T <- T + 1
            }
          }
        }
        tmp1 <- c(tmp1,tmp3)
        tmp2 <- c(tmp2,tmp4)
      }
      pl.res[(cc[s]+1):cc[s+1],1] <- tmp1
      pl.res[(cc[s]+1):cc[s+1],2] <- tmp2
    }
  }
  if(Rem == 1){
    yr1 <- rbinom(1,1,ifelse(Ytil1[N]==1,Se[1],1-Sp[1]))
    yr2 <- rbinom(1,1,ifelse(Ytil2[N]==1,Se[2],1-Sp[2]))
    Zmat <- rbind(Zmat,c(yr1,yr2,1,Se[1],Se[2],
                 Sp[1],Sp[2],N,rep(-9,c.s[1]-1)))
    T <- T + 1
  }
  if(Rem > 1){
    rid <- (M*c.s[1]+1):N
    ytr1 <- Ytil1[rid]
    ytr2 <- Ytil2[rid]
    zr1 <- rbinom(1,1,ifelse(sum(ytr1)>0,Se[1],1-Sp[1]))
    zr2 <- rbinom(1,1,ifelse(sum(ytr2)>0,Se[2],1-Sp[2]))
    Zmat <- rbind(Zmat,c(zr1,zr2,Rem,Se[1],Se[2],
                  Sp[1],Sp[2],rid,rep(-9,c.s[1]-Rem)))
    T <- T + 1
    if(zr1+zr2 > 0){
      yrm1 <- rbinom(Rem,1,ifelse(ytr1==1,Se[1],1-Sp[1]))
      yrm2 <- rbinom(Rem,1,ifelse(ytr2==1,Se[2],1-Sp[2]))
	  Zmat <- rbind(Zmat,cbind(yrm1,yrm2,1,Se[1],Se[2],Sp[1],
	                Sp[2],rid,matrix(-9,Rem,c.s[1]-1)))
      T <- T + Rem
    }
  }
  ivid <- paste( rep("Indv",psz[[1]][1]),1:psz[[1]][1],sep="" )
  colnames(Zmat) <- c("Z1","Z2","psz","Se1","Se2","Sp1","Sp2",ivid)
  rownames(Zmat) <- paste("Pool:",1:nrow(Zmat),sep="")
  return(Zmat)
}

#******************************************************#
# Support function for hier.alg.data()
#
mpt.data <- function(p,N,design,Se,Sp,Yt){
  J <- 2
  psz <- design
  
  id <- sample( 1:N, N, replace=FALSE )
  Yt <- Yt[ id, ]

  if(!is.null(Yt)){
    Ytil1 <- Yt[ ,1]
    Ytil2 <- Yt[ ,2]
  }
  else{
    Ymul <- rmultinom(N,1,p)
    Ytil1 <- Ymul[2, ] + Ymul[4, ]
    Ytil2 <- Ymul[3, ] + Ymul[4, ]
  }  

  tmp <- floor(N/psz)
  Rem <- N-psz*tmp
  cvec <- rep(psz, tmp)
  if( Rem > 0 ){
    cvec <- c( rep(psz, tmp), Rem )
  }
  M <- length( cvec )
  Zmat <- matrix(-9, M, 7+psz)
  id <- cumsum(c(0,cvec))
  pl.res <- matrix(-9,M,2)
  for(m in 1:M){
    mid <- (id[m]+1):id[m+1]
    prob1 <- ifelse(sum(Ytil1[mid])>0,Se[1],1-Sp[1])
    prob2 <- ifelse(sum(Ytil2[mid])>0,Se[2],1-Sp[2])
    Zmat[m,1] <- rbinom(1,1,prob1)
    Zmat[m,2] <- rbinom(1,1,prob2)
    Zmat[m,3] <- length( mid )
    Zmat[m,4:7] <- c(Se,Sp)
    Zmat[m,8:(7+length(mid))] <- mid
  }

  ivid <- paste( rep("Indv",psz),1:psz,sep="" )
  colnames(Zmat) <- c("Z1","Z2","psz","Se1","Se2","Sp1","Sp2",ivid)
  rownames(Zmat) <- paste("Pool:",1:nrow(Zmat),sep="")
  return(Zmat)
}

#******************************************************#
# Main function for simulating 2-dim array pooling data with L = 1 assay
#
array.2dim.data <- function(N,p=NULL,design,Se,Sp,Yt=NULL){
  n <- design[1]
  L <- floor(N/n^2)
  N0 <- L*n^2
  Rem <- N-N0
  npools <- L*2*n
  
  id <- sample( 1:N, N, replace=FALSE )
  Yt <- Yt[ id, ]

  if(!is.null(Yt)){
    Ytil1 <- Yt[ ,1]
    Ytil2 <- Yt[ ,2]
  }
  else{
    Ymul <- rmultinom(N,1,p)
    Ytil1 <- Ymul[2, ] + Ymul[4, ]
    Ytil2 <- Ymul[3, ] + Ymul[4, ]
  } 

  id1 <- cumsum( c(0,rep(n^2,L)) )
  Zmat <- Ymat <- NULL
  tmp.id <- rep(-9,N0)
  ntest <- count2 <- 0
  for(l in 1:L){
    Z_id <- matrix((id1[l]+1):id1[l+1],n,n)
    tmp1 <- Ytil1[ (id1[l]+1):id1[l+1] ]
    Ymat1 <- matrix(tmp1,n,n)
    tmp2 <- Ytil2[ (id1[l]+1):id1[l+1] ]
    Ymat2 <- matrix(tmp2,n,n)
    R1 <- rbinom(n,1,ifelse(rowSums(Ymat1)>0,Se[1],1-Sp[1]))
    R2 <- rbinom(n,1,ifelse(rowSums(Ymat2)>0,Se[2],1-Sp[2]))
    Zmat <- rbind(Zmat,cbind(R1,R2,n,Se[1],Se[2],Sp[1],Sp[2],Z_id))
    C1 <- rbinom(n,1,ifelse(colSums(Ymat1)>0,Se[1],1-Sp[1]))
    C2 <- rbinom(n,1,ifelse(colSums(Ymat2)>0,Se[2],1-Sp[2]))
    Zmat <- rbind(Zmat,cbind(C1,C2,n,Se[1],Se[2],Sp[1],Sp[2],t(Z_id)))
    ntest <- ntest + 2*n
    count <- 0
    tpid <- matrix(-9,n,n)
    for(i in 1:n){
      for(j in 1:n){
	    T1 <- 0
	    if(R1[i]==1 & C1[j]==1) T1 <- T1 + 1
	    if(R1[i]==1 & sum(C1)==0) T1 <- T1 + 1
	    if(sum(R1)==0 & C1[j]==1) T1 <- T1 + 1
	    T2 <- 0
	    if(R2[i]==1 & C2[j]==1) T2 <- T2 + 1
	    if(R2[i]==1 & sum(C2)==0) T2 <- T2 + 1
	    if(sum(R2)==0 & C2[j]==1) T2 <- T2 + 1
	    if(T1+T2>=1){
          prob1 <- ifelse( Ymat1[i,j]>0, Se[1], 1-Sp[1] )
          y1 <- rbinom(1,1,prob1)
          prob2 <- ifelse( Ymat2[i,j]>0, Se[2], 1-Sp[2] )
          y2 <- rbinom(1,1,prob2)
		  count <- count + 1
		  Ymat <- rbind(Ymat,c(y1,y2,1,Se[1],Se[2],Sp[1],Sp[2],Z_id[i,j]))
		  count2 <- count2 + 1
		  tpid[i,j] <- npools + count2
	    }
	  }
    }
    tmp.id[ (id1[l]+1):id1[l+1] ] <- tpid
    ntest <- ntest + count
  }
  Zmat <- rbind(Zmat,cbind(Ymat,matrix(-9,nrow(Ymat),n-1)))
  plrem <- Zrem <- NULL
  if(Rem >= 1){
	if(Rem > n){
	  L0 <- floor(Rem/n)
	  N1 <- Rem-L0*n
	  if(N1 > 0){
	    id2 <- cumsum(c(0,rep(n,L0),N1))
	  }else{
	    id2 <- cumsum(c(0,rep(n,L0)))
	  }
	  for(r in 1:(length(id2)-1)){
	     Id1 <- N0+id2[r]+1
		 Id2 <- N0+id2[r+1]
		 yrem1 <- Ytil1[Id1:Id2]
		 yrem2 <- Ytil2[Id1:Id2]
	     res <- array.rem.test(yrem1,yrem2,Se,Sp,Id1,Id2,n,ntest)
	     Zrem <- rbind(Zrem,res$Zrem)
		 plrem <- rbind(plrem,res$plrem)
		 ntest <- ntest + res$trem
	  }
    }else{
	  Id1 <- N0+1
	  Id2 <- N
	  yrem1 <- Ytil1[Id1:Id2]
	  yrem2 <- Ytil2[Id1:Id2]
	  res <- array.rem.test(yrem1,yrem2,Se,Sp,Id1,Id2,n,ntest)
	  Zrem <- rbind(Zrem,res$Zrem)
	  plrem <- rbind(plrem,res$plrem)
	  ntest <- ntest + res$trem
	}
  }
  Zmat <- rbind(Zmat,Zrem)
  rownames(Zmat) <- NULL
  idv.id <- paste( rep("Indv",n), 1:n, sep="" )
  colnames(Zmat) <- c("Z1","Z2","psz","Se1","Se2","Sp1","Sp2",idv.id)
  rownames(Zmat) <- paste("Pool:",1:nrow(Zmat),sep="")
  T <- nrow(Zmat)
  return(Zmat)
}

#******************************************************#
# Support function for array.2dim.data()
#
array.rem.test <- function(yrem1,yrem2,Se,Sp,Id1,Id2,n,T){
  psz <- length(yrem1)
  if(psz == 1){
	 yr1 <- rbinom(1,1,ifelse(yrem1==1,Se[1],1-Sp[1]))
	 yr2 <- rbinom(1,1,ifelse(yrem2==1,Se[2],1-Sp[2]))
     Zrem <- c(yr1,yr2,1,Se[1],Se[2],Sp[1],Sp[2],Id1,rep(-9,n-1))
	 count3 <- 1
	 plrem <- matrix(c(T+1,rep(-9,2)),ncol=3)
  }
  if(psz > 1){
    Yrem <- NULL
    psz.id <- Id1:Id2
    zrem1 <- rbinom(1,1,ifelse(sum(yrem1)>0,Se[1],1-Sp[1]))
    zrem2 <- rbinom(1,1,ifelse(sum(yrem2)>0,Se[2],1-Sp[2]))
    Zrem <- c(zrem1,zrem2,psz,Se[1],Se[2],Sp[1],Sp[2],psz.id)
    count3 <- 1
    fill1 <- matrix(-9,psz,2)
    if(zrem1+zrem2 > 0){
	  yrem1 <- rbinom(psz,1,ifelse(yrem1>0,Se[1],1-Sp[1]))
	  yrem2 <- rbinom(psz,1,ifelse(yrem2>0,Se[2],1-Sp[2]))
      Yrem <- cbind(yrem1,yrem2,1,Se[1],Se[2],Sp[1],Sp[2],psz.id,matrix(-9,psz,psz-1))
	  count3 <- count3 + psz
	  fill1 <- cbind(T+1+(1:psz),rep(-9,psz))
    }
	plrem <- cbind(rep(T+1,psz),fill1)
	Zrem <- rbind(Zrem,Yrem)
	if(psz < n) Zrem <- cbind(Zrem,matrix(-9,nrow(Zrem),n-psz))
  }
  list("Zrem"=Zrem,"plrem"=plrem,"trem"=count3)
}

#******************************************************#
# Main function for simulating hierarchical pooling data with L = 2 assays
#
hier.alg.data_L2 <- function(N,p=NULL,design,Se,Sp,Yt=NULL){
  J <- 2
  c.s <- design
  S <- length(c.s)
  
  id <- sample( 1:N, N, replace=FALSE )
  Yt <- Yt[ id, ]

  if(!is.null(Yt)){
    Ytil1 <- Yt[ ,1]
    Ytil2 <- Yt[ ,2]
  }
  else{
    Ymul <- rmultinom(N,1,p)
    Ytil1 <- Ymul[2, ] + Ymul[4, ]
    Ytil2 <- Ymul[3, ] + Ymul[4, ]
  }  

  Ytmat <- matrix(-9,N,3+S)
  Ytmat[ ,1] <- Ytil1
  Ytmat[ ,2] <- Ytil2
  M <- floor(N/c.s[1])
  Rem <- N-M*c.s[1]
  psz <- n.div <- list()
  psz[[1]] <- rep(c.s[1],M)
  n.div[[1]] <- rep(1,length(psz[[1]]))
  n.sm <- matrix(-9,length(psz[[1]]),S)
  n.sm[ ,1] <- 1
  if(S > 1){
    for( s in 1:(S-1) ){
	  store <- tmp <- NULL
	  for(i in 1:length(psz[[s]])){
	    temp <- rep( c.s[s+1], floor(psz[[s]][i]/c.s[s+1]) )
  	    store <- c(store,temp)
	    tmp <- c(tmp,length(temp))
	  }
	  psz[[s+1]] <- store
	  n.div[[s+1]] <- tmp
    }
    vec <- rep(1,length(psz[[1]]))
    for(s in 1:(S-1) ){
	  id0 <- cumsum(c(0,vec))
	  for(i in 1:length(psz[[1]])){
        vec[i] <- sum(n.div[[s+1]][(id0[i]+1):id0[i+1]])
      }
      n.sm[ ,s+1] <- vec
    }
  }
  T <- 0
  Zmat <- NULL
  cc <- cumsum(c(0,colSums(n.sm)))
  id <- cumsum(c(0,psz[[1]]))
  pl.res <- matrix(-9,cc[S+1],2)
  for(m in 1:M){
    mid <- (id[m]+1):id[m+1]
    prob1 <- ifelse(sum(Ytil1[mid])>0,Se[1,1],1-Sp[1,1])
    prob2 <- ifelse(sum(Ytil2[mid])>0,Se[1,2],1-Sp[1,2])
    z1 <- rbinom(1,1,prob1)
    z2 <- rbinom(1,1,prob2)
    pl.res[m,1] <- z1
    pl.res[m,2] <- z2
    Zmat <- rbind(Zmat,c(z1,z2,length(mid),Se[1,1],Se[1,2],Sp[1,1],Sp[1,2],mid))
    Ytmat[mid,J+1+1] <- m
    T <- T + 1
  }
  if( S > 1){
    for(s in 2:S){
      Z1 <- pl.res[(cc[s-1]+1):cc[s],1]
      Z2 <- pl.res[(cc[s-1]+1):cc[s],2]
      cid <- cumsum(c(0,psz[[s]]))
      cn <- cumsum(c(0,n.div[[s]]))
      tmp1 <- tmp2 <- NULL
      for(d in 1:length(psz[[s-1]])){
        tmp3 <- tmp4 <- NULL
        if(Z1[d]==0 & Z2[d]==0){
          tmp3 <- tmp4 <- rep(0,length((cn[d]+1):cn[d+1]))
        }
        if(psz[[s-1]][d]==1){
          tmp3 <- Z1[d]
          tmp4 <- Z2[d]
        }
        if(psz[[s-1]][d]>1){
          if(Z1[d]+Z2[d] > 0){
            for(i in (cn[d]+1):cn[d+1]){
              crng <- (cid[i]+1):cid[i+1]
              prob1 <- ifelse(sum(Ytil1[crng])>0,Se[s,1],1-Sp[s,1])
              prob2 <- ifelse(sum(Ytil2[crng])>0,Se[s,2],1-Sp[s,2])
              ztp1 <- rbinom(1,1,prob1)
              ztp2 <- rbinom(1,1,prob2)
              tmp3 <- c(tmp3,ztp1)
              tmp4 <- c(tmp4,ztp2)
              fill1 <- rep(-9,psz[[1]][1]-length(crng))
              Zmat <- rbind(Zmat,c(ztp1,ztp2,length(crng),Se[s,1],
                            Se[s,2],Sp[s,1],Sp[s,2],crng,fill1))
              T <- T + 1
              Ytmat[crng,4+s-1] <- T
            }
          }
        }
        tmp1 <- c(tmp1,tmp3)
        tmp2 <- c(tmp2,tmp4)
      }
      pl.res[(cc[s]+1):cc[s+1],1] <- tmp1
      pl.res[(cc[s]+1):cc[s+1],2] <- tmp2
    }
  }
  if(Rem == 1){
    yr1 <- rbinom(1,1,ifelse(Ytil1[N]==1,Se[S,1],1-Sp[S,1]))
    yr2 <- rbinom(1,1,ifelse(Ytil2[N]==1,Se[S,2],1-Sp[S,2]))
    Zmat <- rbind(Zmat,c(yr1,yr2,1,Se[S,1],Se[S,2],
                 Sp[S,1],Sp[S,2],N,rep(-9,c.s[1]-1)))
    T <- T + 1
    Ytmat[N,J+1+1] <- T
  }
  if(Rem > 1){
    rid <- (M*c.s[1]+1):N
    ytr1 <- Ytil1[rid]
    ytr2 <- Ytil2[rid]
    zr1 <- rbinom(1,1,ifelse(sum(ytr1)>0,Se[S-1,1],1-Sp[S-1,1]))
    zr2 <- rbinom(1,1,ifelse(sum(ytr2)>0,Se[S-1,2],1-Sp[S-1,2]))
    Zmat <- rbind(Zmat,c(zr1,zr2,Rem,Se[S-1,1],Se[S-1,2],
                  Sp[S-1,1],Sp[S-1,2],rid,rep(-9,c.s[1]-Rem)))
    T <- T + 1
    Ytmat[rid,J+1+1] <- T
    if(zr1+zr2 > 0){
      yrm1 <- rbinom(Rem,1,ifelse(ytr1==1,Se[S,1],1-Sp[S,1]))
      yrm2 <- rbinom(Rem,1,ifelse(ytr2==1,Se[S,2],1-Sp[S,2]))
	  Zmat <- rbind(Zmat,cbind(yrm1,yrm2,1,Se[S,1],Se[S,2],Sp[S,1],
	                Sp[S,2],rid,matrix(-9,Rem,c.s[1]-1)))
      Ytmat[rid,J+1+2] <- (T+1):(T+Rem)
      T <- T + Rem
    }
  }
  temp <- Ytmat[ ,-(1:3)]
  for(j in 1:N) Ytmat[j,3] <- sum(temp[j, ]>0)
  ivid <- paste( rep("Indv",psz[[1]][1]),1:psz[[1]][1],sep="" )
  
  uniques <- unique(Zmat[ ,4])
  L <- length(uniques)
  Lid <- rep(-9,nrow(Zmat))
  for(l in 1:L){
	Mid <- uniques[l]==Zmat[ ,4]
	Lid[Mid] <- l
  }
  Zmat2 <- cbind(Zmat[ ,1:7], Lid, Zmat[ ,-(1:7)])

  colnames(Zmat) <-  c("Z1","Z2","psz","Se1","Se2","Sp1","Sp2",ivid)
  colnames(Zmat2) <- c("Z1","Z2","psz","Se1","Se2","Sp1","Sp2","Assay",ivid)  
#  rownames(Zmat) <- paste("Pool:",1:nrow(Zmat),sep="")
  T <- nrow(Zmat2)
  return( Zmat2 )
}

#******************************************************#
# Main function for simulating 2-dim array pooling data with L = 2 assays
#
array.2dim.data_L2 <- function(N,p=NULL,design,Se,Sp,Yt=NULL){
  n <- design[1]
  L <- floor(N/n^2)
  N0 <- L*n^2
  Rem <- N-N0
  npools <- L*2*n

  id <- sample( 1:N, N, replace=FALSE )
  Yt <- Yt[ id, ]

  if(!is.null(Yt)){
    Ytil1 <- Yt[ ,1]
    Ytil2 <- Yt[ ,2]
  }
  else{
    Ymul <- rmultinom(N,1,p)
    Ytil1 <- Ymul[2, ] + Ymul[4, ]
    Ytil2 <- Ymul[3, ] + Ymul[4, ]
  }   

  id1 <- cumsum( c(0,rep(n^2,L)) )
  Zmat <- Ymat <- NULL
  tmp.id <- rep(-9,N0)
  ntest <- count2 <- 0
  for(l in 1:L){
    Z_id <- matrix((id1[l]+1):id1[l+1],n,n)
    tmp1 <- Ytil1[ (id1[l]+1):id1[l+1] ]
    Ymat1 <- matrix(tmp1,n,n)
    tmp2 <- Ytil2[ (id1[l]+1):id1[l+1] ]
    Ymat2 <- matrix(tmp2,n,n)
    R1 <- rbinom(n,1,ifelse(rowSums(Ymat1)>0,Se[1,1],1-Sp[1,1]))
    R2 <- rbinom(n,1,ifelse(rowSums(Ymat2)>0,Se[1,2],1-Sp[1,2]))
    Zmat <- rbind(Zmat,cbind(R1,R2,n,Se[1,1],Se[1,2],Sp[1,1],Sp[1,2],Z_id))
    C1 <- rbinom(n,1,ifelse(colSums(Ymat1)>0,Se[1,1],1-Sp[1,1]))
    C2 <- rbinom(n,1,ifelse(colSums(Ymat2)>0,Se[1,2],1-Sp[1,2]))
    Zmat <- rbind(Zmat,cbind(C1,C2,n,Se[1,1],Se[1,2],Sp[1,1],Sp[1,2],t(Z_id)))
    ntest <- ntest + 2*n
    count <- 0
    tpid <- matrix(-9,n,n)
    for(i in 1:n){
      for(j in 1:n){
	    T1 <- 0
	    if(R1[i]==1 & C1[j]==1) T1 <- T1 + 1
	    if(R1[i]==1 & sum(C1)==0) T1 <- T1 + 1
	    if(sum(R1)==0 & C1[j]==1) T1 <- T1 + 1
	    T2 <- 0
	    if(R2[i]==1 & C2[j]==1) T2 <- T2 + 1
	    if(R2[i]==1 & sum(C2)==0) T2 <- T2 + 1
	    if(sum(R2)==0 & C2[j]==1) T2 <- T2 + 1
	    if(T1+T2>=1){
          prob1 <- ifelse( Ymat1[i,j]>0, Se[2,1], 1-Sp[2,1] )
          y1 <- rbinom(1,1,prob1)
          prob2 <- ifelse( Ymat2[i,j]>0, Se[2,2], 1-Sp[2,2] )
          y2 <- rbinom(1,1,prob2)
		  count <- count + 1
		  Ymat <- rbind(Ymat,c(y1,y2,1,Se[2,1],Se[2,2],Sp[2,1],Sp[2,2],Z_id[i,j]))
		  count2 <- count2 + 1
		  tpid[i,j] <- npools + count2
	    }
	  }
    }
    tmp.id[ (id1[l]+1):id1[l+1] ] <- tpid
    ntest <- ntest + count
  }
  Zmat <- rbind(Zmat,cbind(Ymat,matrix(-9,nrow(Ymat),n-1)))
  plrem <- Zrem <- NULL
  if(Rem >= 1){
	if(Rem > n){
	  L0 <- floor(Rem/n)
	  N1 <- Rem-L0*n
	  if(N1 > 0){
	    id2 <- cumsum(c(0,rep(n,L0),N1))
	  }else{
	    id2 <- cumsum(c(0,rep(n,L0)))
	  }
	  for(r in 1:(length(id2)-1)){
	     Id1 <- N0+id2[r]+1
	     Id2 <- N0+id2[r+1]
	     yrem1 <- Ytil1[Id1:Id2]
	     yrem2 <- Ytil2[Id1:Id2]
	     res <- array.rem.test(yrem1,yrem2,Se,Sp,Id1,Id2,n,ntest)
	     Zrem <- rbind(Zrem,res$Zrem)
	     plrem <- rbind(plrem,res$plrem)
	     ntest <- ntest + res$trem
	  }
    }else{
	  Id1 <- N0+1
	  Id2 <- N
	  yrem1 <- Ytil1[Id1:Id2]
	  yrem2 <- Ytil2[Id1:Id2]
	  res <- array.rem.test(yrem1,yrem2,Se,Sp,Id1,Id2,n,ntest)
	  Zrem <- rbind(Zrem,res$Zrem)
	  plrem <- rbind(plrem,res$plrem)
	  ntest <- ntest + res$trem
	}
  }
  Zmat <- rbind(Zmat,Zrem)
  rownames(Zmat) <- NULL
  ivid <- paste( rep("Indv",n), 1:n, sep="" )

  uniques <- unique(Zmat[ ,4])
  L <- length(uniques)
  Lid <- rep(-9,nrow(Zmat))
  for(l in 1:L){
	Mid <- uniques[l]==Zmat[ ,4]
	Lid[Mid] <- l
  }
  Zmat2 <- cbind(Zmat[ ,1:7], Lid, Zmat[ ,-(1:7)])

  colnames(Zmat) <-  c("Z1","Z2","psz","Se1","Se2","Sp1","Sp2",ivid)
  colnames(Zmat2) <- c("Z1","Z2","psz","Se1","Se2","Sp1","Sp2","Assay",ivid)  
#  rownames(Zmat) <- paste("Pool:",1:nrow(Zmat),sep="")
  T <- nrow(Zmat2)
  return(Zmat2)
}

#******************************************************#
# Support function for array.2dim.data_L2()
#
array.rem.test_L2 <- function(yrem1,yrem2,Se,Sp,Id1,Id2,n,T){
  psz <- length(yrem1)
  if(psz == 1){
	 yr1 <- rbinom(1,1,ifelse(yrem1==1,Se[2,1],1-Sp[2,1]))
	 yr2 <- rbinom(1,1,ifelse(yrem2==1,Se[2,2],1-Sp[2,2]))
       Zrem <- c(yr1,yr2,1,Se[2,1],Se[2,2],Sp[2,1],Sp[2,2],Id1,rep(-9,n-1))
	 count3 <- 1
	 plrem <- matrix(c(T+1,rep(-9,2)),ncol=3)
  }
  if(psz > 1){
    Yrem <- NULL
    psz.id <- Id1:Id2
    zrem1 <- rbinom(1,1,ifelse(sum(yrem1)>0,Se[1,1],1-Sp[1,1]))
    zrem2 <- rbinom(1,1,ifelse(sum(yrem2)>0,Se[1,2],1-Sp[1,2]))
    Zrem <- c(zrem1,zrem2,psz,Se[1,1],Se[1,2],Sp[1,1],Sp[1,2],psz.id)
    count3 <- 1
    fill1 <- matrix(-9,psz,2)
    if(zrem1+zrem2 > 0){
	  yrem1 <- rbinom(psz,1,ifelse(yrem1>0,Se[2,1],1-Sp[2,1]))
	  yrem2 <- rbinom(psz,1,ifelse(yrem2>0,Se[2,2],1-Sp[2,2]))
        Yrem <- cbind(yrem1,yrem2,1,Se[2,1],Se[2,2],Sp[2,1],Sp[2,2],psz.id,matrix(-9,psz,psz-1))
	  count3 <- count3 + psz
	  fill1 <- cbind(T+1+(1:psz),rep(-9,psz))
    }
    plrem <- cbind(rep(T+1,psz),fill1)
    Zrem <- rbind(Zrem,Yrem)
    if(psz < n) Zrem <- cbind(Zrem,matrix(-9,nrow(Zrem),n-psz))
  }
  list("Zrem"=Zrem,"plrem"=plrem,"trem"=count3)
}

