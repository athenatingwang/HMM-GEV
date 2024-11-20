##Function for ordinary pseudo-residuals
HMMNLM.resid <- function(data,HMMest,distn)
{
  nn <- length(data)
  m<-length(HMMest$delta)
  gammah <- HMMest$Pi
  deltah <- HMMest$delta
  if(distn=="norm"){
    muh <- HMMest$mu
    sigh <- HMMest$sigma
    pRS<-matrix(as.double(0), nrow = nn, ncol = m)
    for (k in 1:m) pRS[, k] <- dnorm(x = data,mean=muh[k],sd=sigh[k])
  }
  if(distn=="lnorm"){
    muh <- HMMest$mu
    sigh <- HMMest$sigma
    pRS<-matrix(as.double(0), nrow = nn, ncol = m)
    for (k in 1:m) pRS[, k] <- dlnorm(x = data,meanlog=muh[k],sdlog=sigh[k])
  }
  if(distn=="weibull"){
    shapeh <- HMMest$shape
    rateh <- HMMest$rate
    pRS<-matrix(as.double(0), nrow = nn, ncol = m)
    for (k in 1:m) pRS[, k] <- dweibull(x= data,shape=shapeh[k],scale=1/rateh[k])
  }
  if(distn=="gamma"){
    shapeh <- HMMest$shape
    rateh <- HMMest$rate
    pRS<-matrix(as.double(0), nrow = nn, ncol = m)
    for (k in 1:m) pRS[, k] <- dgamma(x = data,shape=shapeh[k],rate=rateh[k])
  }
  if(distn=="gev"){
    loch <- HMMest$loct
    scaleh <- HMMest$scale
    shapeh <- HMMest$shape
    pRS<-matrix(as.double(0), nrow = nn, ncol = m)
    for (k in 1:m) pRS[,k] <- devd(x = data,loc=loch[k], scale=scaleh[k], 
                                   shape=shapeh[k],log=FALSE,type = "GEV")
  }
  
  
  #scaled forward variable
  logalpha <- matrix(as.double(rep(0, nn * m)), nrow = nn)
  lscale <- as.double(0)
  phi <- as.double(deltah)
  for (t in 1:nn)
  {
    if (t > 1)
      phi <- phi %*% gammah
    phi <- phi * pRS[t,]
    sumphi <- sum(phi)
    phi <- phi/sumphi
    lscale <- lscale + log(sumphi)
    logalpha[t, ] <- log(phi) + lscale
  }
  # LL <- lscale
  #
  ##Scaled backward variable
  logbeta <- matrix(as.double(rep(0, nn * m)), nrow = nn)
  phi <- as.double(rep(1/m, m))
  lscale <- as.double(log(m))
  for (t in seq(nn - 1, 1, -1))
  {
    phi <- gammah %*% (pRS[t+1,] * phi)
    logbeta[t, ] <- log(phi) + lscale
    sumphi <- sum(phi)
    phi <- phi/sumphi
    lscale <- lscale + log(sumphi)
  }
  if(distn=="norm"){
    for (i in 1:nn) pRS[i, ] <- pnorm(data[i], mean=muh,sd=sigh)
  }
  if(distn=="lnorm"){
    for (i in 1:nn) pRS[i, ] <- plnorm(data[i], meanlog=muh,sdlog=sigh)
  }
  if(distn=="weibull"){
    for (i in 1:nn) pRS[i, ] <- pweibull(data[i], shape=shapeh,scale=1/rateh)
  }
  if(distn=="gamma"){
    for (i in 1:nn) pRS[i, ] <- pgamma(data[i], shape=shapeh,rate=rateh)
  }
  if(distn=="gev"){
    for (k in 1:m) pRS[,k] <- pevd(data,loc=loch[k], scale=scaleh[k], 
                                   shape=shapeh[k],log=FALSE,type = "GEV")
  }
  prob <- rep(NA, nn)
  for (i in 1:nn) {
    if (i == 1) {
      pre <- deltah
    }
    else {
      la <- logalpha[i - 1, ]
      pre <- exp(la - mean(la[la != -Inf])) %*% gammah
    }
    lb <- logbeta[i, ]
    post <- exp(lb - mean(lb[lb != -Inf]))
    prob[i] <- (pre %*% diag(pRS[i, ]) %*% post)/(pre %*%
                                                    post)
  }
  return(qnorm(prob))
}


