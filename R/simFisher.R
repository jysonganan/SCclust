#kXk matrices yy, ny, t(ny) and nn contain matrix elements [1,1], [2,1], [1,2] and [2,2] of 
#the 2X2 contingency tables for all possible combinations of k variables. For a variable i
#return its Fisher exact p-values with variables j>=i. 
fastFisher<-function(i,yy,ny,nn){
  ftp<-rep(1,nrow(yy))
  for(j in i:nrow(yy))ftp[j]<-fisher.test(matrix(nrow=2,ncol=2,
                                                 data=c(yy[i,j],ny[i,j],ny[j,i],nn[i,j])),alternative="greater")$p.value
  return(ftp)
}

#Metropolis randomization of vector x of 0s and 1s, where p is a vector of occupancy
#probabilities, keeping the sum of x fixed. If there are more 0s than 1s in x, a subset of
#vacant positions is chosen at random, and 1s are moved to these positions with Metropolis
#transition probability. If there are more 1s than 0s, the roles of 0s and 1s are inetrchanged
metro<-function(x,p,sweeps){
  for(i in 1:sweeps){
    occ<-which(x==1&p<1)
    vac<-which(x==0&p>0)
    locc<-length(occ)
    lvac<-length(vac)
    if(locc==0|lvac==0)return(x)
    if(locc>lvac){
      orig<-sample(occ,size=lvac)
      dest<-vac
    }
    else{
      orig<-occ
      dest<-sample(vac,size=locc)
    }
    lmove<-length(orig)
    moves<-(p[dest]*(1-p[orig])/p[orig]/(1-p[dest])>runif(length(orig)))
    x[dest[moves]]<-1
    x[orig[moves]]<-0
  }
  return(x)
}

#Stouffer and Fisher combos come courtesy of Wikipedia: 
#http://en.wikipedia.org/wiki/Fisher's_method . I only changed the function names.
StoufferCombo <- function(p, w) { # p is a vector of p-values
  if (missing(w)) {
    w <- rep(1, length(p))/length(p)
  } else {
    if (length(w) != length(p))
      stop("Length of p and w must equal!")
  }
  Zi <- qnorm(1-p) 
  Z  <- sum(w*Zi)/sqrt(sum(w^2))
  p.val <- 1-pnorm(Z)
  return(c(Z = Z, p.value = p.val))
}

FisherCombo <- function(p) {
  Xsq <- -2*sum(log(p))
  p.val <- pchisq(Xsq, df = 2*length(p),lower.tail=FALSE)
  return(c(Xsq = Xsq, p.value = p.val))
}

simFisher<-function(m,nsim,nsweep,seedme,distrib=c("vanilla","Rparallel"),njobs=1,
                    combo=c("Fisher","Stouffer")){
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seedme)
  distrib<-match.arg(distrib)
  combo<-match.arg(combo)
  if(distrib=="Rparallel"){
    ncores<-min(njobs,parallel::detectCores())
    cl<-parallel::makeCluster(getOption("cl.cores",ncores))
    parallel::clusterSetRNGStream(cl)
  }
  if(!is.list(m))m<-list(m)
  rf<-lapply(m,rowMeans)
  tp<-matrix(ncol=nsim,nrow=ncol(m[[1]])*(ncol(m[[1]])-1)/2)
  xmat<-matrix(ncol=length(m),nrow=ncol(m[[1]])*(ncol(m[[1]])-1)/2)
  for(i in 1:nsim){
    for(j in 1:length(m)){
      if(nsweep>0){
        if(distrib=="vanilla")m[[j]]<-apply(m[[j]],2,metro,p=rf[[j]],sweeps=nsweep)
        if(distrib=="Rparallel")m[[j]]<-
            parallel::parApply(cl=cl,X=m[[j]],MARGIN=2,FUN=metro,p=rf[[j]],sweeps=nsweep)
      }
      yy<-t(m[[j]])%*%m[[j]]
      ny<-t(1-m[[j]])%*%m[[j]]
      nn<-t(1-m[[j]])%*%(1-m[[j]])
      if(distrib=="vanilla")x[[j]]<-sapply(1:nrow(yy),FUN=fastFisher,yy=yy,ny=ny,nn=nn)
      if(distrib=="Rparallel"){
        lbi<-1:nrow(yy)
        lbi[lbi%%2==0]<-nrow(yy)+2-nrow(yy)%%2-lbi[lbi%%2==0]
        x<-parallel::parSapply(cl,X=lbi,FUN=fastFisher,yy=yy,ny=ny,nn=nn)[,order(lbi)]
      }
      x<-pmin(x,t(x))
      xmat[,j]<-x[upper.tri(x)]	
    }
    if(length(m)==1)tp[,i]<-xmat[,1]
    else{ 
      if(combo=="Fisher"){
        Xsq <- -2*rowSums(log(xmat))
        tp[,i]<-sapply(Xsq,pchisq,df=2*ncol(xmat),lower.tail=FALSE)
      }
      if(combo=="Stouffer"){
        Z<-colSums(apply(1-xmat,1,qnorm))/sqrt(ncol(xmat))
        tp[,i]<-1-sapply(Z,pnorm)
      }
    }
  }	
  parallel::stopCluster(cl)
  #dump results
  return(tp)
}
