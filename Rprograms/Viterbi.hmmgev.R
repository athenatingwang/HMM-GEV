Viterbi.hmmgev <- function(data, HMMest){
  n <- length(data)
  loch <- HMMest$loct
  scaleh <- HMMest$scale
  shapeh <- HMMest$shape
  gamman <- HMMest$Pi
  delta <- HMMest$delta
  m <- length(shapeh)
  
  pRS<-matrix(0,n,m)
  pRS <- matrix(as.double(pRS),nrow=n)  
    for(i in 1:m){
      pRS[,i] <- devd(x=data,loc=loch[i], scale=scaleh[i], 
                      shape=shapeh[i],log=FALSE,type = "GEV")
    }

    nu <- matrix(NA, nrow = n, ncol = m)
    y <- rep(NA, n)
    nu[1, ] <- log(delta) + log(pRS[1,])
    logPi <- log(gamman)
    for (i in 2:n) {
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        nu[i, ] <- apply(matrixnu + logPi, 2, max) + log(pRS[i,])
    }
#    if (any(nu[n, ] == -Inf)) 
 #       stop("Problems With Underflow")
    y[n] <- which.max(nu[n, ])
    for (i in seq(n - 1, 1, -1)) y[i] <- which.max(logPi[, y[i + 
        1]] + nu[i, ])
    logalpha <- matrix(as.double(rep(0, n * m)), nrow = n)
    lscale <- as.double(0)
    phi <- as.double(delta)
    for (t in 1:n)
    {
      if (t > 1) 
         phi <- phi %*% gamman
      phi <- phi * pRS[t,]
      sumphi <- sum(phi)
      phi <- phi/sumphi
      lscale <- lscale + log(sumphi)
      logalpha[t, ] <- log(phi) + lscale
    }
    LL <- lscale
#
##Scaled backward variable
    logbeta <- matrix(as.double(rep(0, n * m)), nrow = n)
    phi <- as.double(rep(1/m, m))
    lscale <- as.double(log(m))
    for (t in seq(n - 1, 1, -1)) 
    {
       phi <- gamman %*% (pRS[t+1,] * phi)
       logbeta[t, ] <- log(phi) + lscale
       sumphi <- sum(phi)
       phi <- phi/sumphi
       lscale <- lscale + log(sumphi)
    }
#
##E-step
###Calculate v_t(j)
    v <- exp(logalpha + logbeta - LL)
## y is the estimated hidden states 
## v is the estimated probability of each time point being in each state
  return(list(y=y,v=v))
}


