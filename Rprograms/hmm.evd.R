library(Rcpp)
#sourceCpp("../RprogramHMM/prsloop.evd.cpp")
sourceCpp("../RprogramHMM/loop1.cpp")

library(extRemes)

hmm.evd <- function(param,data,nk,cpp=TRUE){
#hmm.evd <- function(param,data,nk){
    n<-length(data)
  locn <- param[1:nk]
  scalen<-exp(param[(nk+1):(2*nk)])
  shapen<-param[(2*nk+1):(3*nk)]
  
  tgamman<-exp(param[(3*nk+1):length(param)])
  tgamma<-diag(nk)
  
  tgamma[!tgamma]<-tgamman
  gamman<-tgamma/apply(tgamma,1,sum)
  
  delta<-solve(t(diag(nk)-gamman+1),rep(1,nk))
  
  pRS<-matrix(0,n,nk)
  pRS <- matrix(as.double(pRS),nrow=n)  
#  if (cpp!=TRUE){
#s1 <- Sys.time()
    for(i in 1:nk){
      pRS[,i] <- devd(data,loc=locn[i], scale=scalen[i], shape=shapen[i],log=FALSE,type = "GEV")
    }
#s2 <- Sys.time()
#s2-s1
#  } else {
#s3 <- Sys.time()
#  pRS <- prsloop(nk,n,data,locn,shapen,scalen)
#s4 <- Sys.time()
#s4-s3
#  }    
  
  logalpha <- matrix(as.double(rep(0, n * nk)), nrow = n)
  lscale <- as.double(0)
  phi <- as.double(delta)
  gamma <- matrix(as.double(gamman),nrow=nk)
  pRS <- matrix(as.double(pRS),nrow=n)
  if (cpp!=TRUE){
    #  loop1 using R code
#s1 <- Sys.time()
    for (t in 1:n)
    {
      if (t > 1) 
        phi <- phi %*% gamma
      phi <- phi * pRS[t,]
      sumphi <- sum(phi)
      phi <- phi/sumphi
      lscale <- lscale + log(sumphi)
      logalpha[t, ] <- log(phi) + lscale
    }
    LL <- lscale
#s2 <- Sys.time()
#s2-s1
  }else{
    if (!is.double(gamma)) stop("gamma is not double precision")
#s3 <- Sys.time()
    loop1 <- loop1(nk, n, phi, pRS, gamma, logalpha)
    logalpha <- loop1$logalp
    LL <- loop1$lscale
#s4 <- Sys.time()
#s4-s3
  }  
  return(-LL)
}

