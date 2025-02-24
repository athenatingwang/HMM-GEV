library(HiddenMarkov)
library(extRemes)
library(lubridate)


load(file="../MaximaData/R1hrmax.image")
R1hrmax <- R1hrmax
len.data <- length(R1hrmax)



#HMM <- get(load("../HMM-GAMgevRes/HMMgevHr6est780.image"))

## Rerun HMMgev to get a better convergence
#source("hmm.evd.R")
nk <- 6
#Pit <- log(HMM$Pi/diag(HMM$Pi))
#Pi0 <- as.vector(Pit[!diag(nk)])
#loct <- HMM$loct+0.02
#shapet <- HMM$shape+0.02
#scalet <- HMM$scale
#param0<-c(loct,log(scalet),shapet,Pi0)

#y <- nlm(hmm.evd,p=param0,data=R1hrmax,nk=nk,hessian = TRUE,iterlim = 10000,steptol = 1e-5,print.level=2)
#loc.est <- y$estimate[1:(nk)]
#scale.est <- exp(y$estimate[(nk+1):(2*nk)])
#shape.est <- y$estimate[(2*nk+1):(3*nk)]

#tgamman<-exp(y$estimate[(3*nk+1):length(y$estimate)])
#tgamma<-diag(nk)

#tgamma[!tgamma]<-tgamman
#gamman<-y$estimate[(3*nk+1):(length(y$estimate)+nk)]<-tgamma/apply(tgamma,1,sum)

#delta<-solve(t(diag(nk)-gamman+1),rep(1,nk))

#HMMest <- list(loct=loc.est, scale=scale.est, shape=shape.est, Pi=gamman,delta=delta, LL=-y$minimum)
#nLL <- -y$minimum

#eval(parse(text=paste('save(HMMest, file="HMMgevHr',nk,'est780-1.image")',sep='')))

HMMest <- get(load("../HMM-GAMgevRes/HMMgevHr6est780-1.image"))


source("Viterbi.hmmgev.R")
source("HMMNLM.resid.R")
source("qq.CI.R")
source("acf.CI.R")



res.hr <- HMMNLM.resid(R1hrmax,HMMest,distn="gev")

#save(res.hr,file="res.hr.image")

#load("res.hr.image")
#res <- res.hr

#qqhmmHr <- qq(res.hr)

load(file="../HMM-GAMgevRes/envelope hourly maxima.RData")
qqhmmHr <- env.hour


postscript(paste("HMMRhrmax",nk,"states-ACF.eps",sep=""),paper="special",
           width=4.5*2/2,height=9)
par(mfrow=c(2,1), mar=c(4.2, 3, 1, 1))

res <- res.hr
qqhmm <- qqhmmHr
ylims<-range(res,qqhmm$low,qqhmm$upp)
plot(qqhmm$med,sort(res),ylab="",ylim=ylims,xlab="Expected",
     cex.axis=1.2,cex.lab=1.2,pch=16,cex=0.5)
title(ylab="Observed", line=2, cex.lab=1.2)
lines(qqhmm$med,qqhmm$med,lwd=0.6)
lines(qqhmm$med,qqhmm$low,col="red",lty=2,lwd=0.6)
lines(qqhmm$med,qqhmm$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(a) Hourly",bty="n",cex=1.5)

p.out<-length(which((sort(res)>qqhmm$upp)|(sort(res)<qqhmm$low)))/length(res)
if(p.out>0){
  message(paste("observations outside envelope (",round(100*p.out,1),"%)",sep=""))
}
if(p.out==0){ message("no observations outside envelope") }


inc<-1.5
lag.max <- 50
acf.hmm <- acf(res,lag.max=lag.max,plot=F)
acf.hmm.ci <- acf.ci(res)
ylims<-range(acf.hmm.ci,acf.hmm$acf,1)
plot(acf.hmm,main="",ci=0,ylim=ylims,cex.axis=1.2,cex.lab=1.2,lwd=0.6,
     xlab="Lag (hour)",ylab="")
lines(acf.hmm.ci$low,lty=2,col="blue")
lines(acf.hmm.ci$upp,lty=2,col="blue")
title(ylab="ACF", line=2, cex.lab=1.2)
legend("topleft","(d) Hourly",bty="n",cex=1.5)
dev.off()













#####  HMM fitted to daily maxima

load(file="../MaximaData/R1daymax.image")
R1hrmax <- R1daymax
len.data <- length(R1hrmax)



HMM <- get(load("../HMM-GAMgevRes/HMMgevDay6est210.image"))

## Rerun HMMgev to get a better convergence
source("hmm.evd.R")
nk <- 6
Pit <- log(HMM$Pi/diag(HMM$Pi))
Pi0 <- as.vector(Pit[!diag(nk)])
loct <- HMM$loct-0.05
shapet <- HMM$shape+0.05
scalet <- HMM$scale
param0<-c(loct,log(scalet),shapet,Pi0)

y <- nlm(hmm.evd,p=param0,data=R1hrmax,nk=nk,hessian = TRUE,iterlim = 10000,steptol = 1e-6,print.level=2)
loc.est <- y$estimate[1:(nk)]
scale.est <- exp(y$estimate[(nk+1):(2*nk)])
shape.est <- y$estimate[(2*nk+1):(3*nk)]

tgamman<-exp(y$estimate[(3*nk+1):length(y$estimate)])
tgamma<-diag(nk)

tgamma[!tgamma]<-tgamman
gamman<-y$estimate[(3*nk+1):(length(y$estimate)+nk)]<-tgamma/apply(tgamma,1,sum)

delta<-solve(t(diag(nk)-gamman+1),rep(1,nk))

HMMest <- list(loct=loc.est, scale=scale.est, shape=shape.est, Pi=gamman,delta=delta, LL=-y$minimum)
nLL <- -y$minimum





source("Viterbi.hmmgev.R")
source("HMMNLM.resid.R")
source("qq.CI.R")
source("acf.CI.R")



res.day <- HMMNLM.resid(R1hrmax,HMMest,distn="gev")


#qqhmmDay <- qq(res.day)

#save(qqhmmDay,file="qqhmmDay.image")
get(load("../HMM-GAMgevRes/qqhmmDay.image"))



postscript(paste("HMMRdaymax",nk,"states-ACF.eps",sep=""),paper="special",
           width=4.5*2/2,height=9)
par(mfrow=c(2,1), mar=c(4.2, 3, 1, 1))

res <- res.day
qqhmm <- qqhmmDay
ylims<-range(res,qqhmm$low,qqhmm$upp)
plot(qqhmm$med,sort(res),ylab="",ylim=ylims,xlab="Expected",
     cex.axis=1.2,cex.lab=1.2,pch=16,cex=0.5)
title(ylab="Observed", line=2, cex.lab=1.2)
lines(qqhmm$med,qqhmm$med,lwd=0.6)
lines(qqhmm$med,qqhmm$low,col="red",lty=2,lwd=0.6)
lines(qqhmm$med,qqhmm$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(b) Daily",bty="n",cex=1.5)

p.out<-length(which((sort(res)>qqhmm$upp)|(sort(res)<qqhmm$low)))/length(res)
if(p.out>0){
  message(paste("observations outside envelope (",round(100*p.out,1),"%)",sep=""))
}
if(p.out==0){ message("no observations outside envelope") }


inc<-1.5
lag.max <- 50
acf.hmm <- acf(res,lag.max=lag.max,plot=F)
acf.hmm.ci <- acf.ci(res)
ylims<-range(acf.hmm.ci,acf.hmm$acf,1)
plot(acf.hmm,main="",ci=0,ylim=ylims,cex.axis=1.2,cex.lab=1.2,lwd=0.6,
     xlab="Lag (day)",ylab="")
lines(acf.hmm.ci$low,lty=2,col="blue")
lines(acf.hmm.ci$upp,lty=2,col="blue")
title(ylab="ACF", line=2, cex.lab=1.2)
legend("topleft","(e) Daily",bty="n",cex=1.5)
dev.off()












####  HMM fitted to monthly maxima

## Load monthly maximum of R1 1994-2019
load(file="../MaximaData/R1.mthmax.RData")
R1hrmax <- d$rmax
len.data <- length(R1hrmax)



HMM <- get(load("../HMM-GAMgevRes/HMMgev4est876.image"))

## Rerun HMMgev to get a better convergence
source("hmm.evd.R")
nk <- 4
Pit <- log(HMM$Pi/diag(HMM$Pi))
Pi0 <- as.vector(Pit[!diag(nk)])
loct <- HMM$loct+0.02
shapet <- HMM$shape+0.02
scalet <- HMM$scale
param0<-c(loct,log(scalet),shapet,Pi0)

y <- nlm(hmm.evd,p=param0,data=R1hrmax,nk=nk,hessian = TRUE,iterlim = 10000,steptol = 1e-4,print.level=2)
loc.est <- y$estimate[1:(nk)]
scale.est <- exp(y$estimate[(nk+1):(2*nk)])
shape.est <- y$estimate[(2*nk+1):(3*nk)]

tgamman<-exp(y$estimate[(3*nk+1):length(y$estimate)])
tgamma<-diag(nk)

tgamma[!tgamma]<-tgamman
gamman<-y$estimate[(3*nk+1):(length(y$estimate)+nk)]<-tgamma/apply(tgamma,1,sum)

delta<-solve(t(diag(nk)-gamman+1),rep(1,nk))

HMMest <- list(loct=loc.est, scale=scale.est, shape=shape.est, Pi=gamman,delta=delta, LL=-y$minimum)
nLL <- -y$minimum



library(mnormt)
gev.param <- y$estimate[1:(3*nk)]
gev.varc <- solve(y$hessian[1:(3*nk),1:(3*nk)])



# Estimated return levels
p50 <- 1/(50*12)
p100 <- 1/(100*12)
p200 <- 1/(200*12)
p500 <- 1/(500*12)


rl.50 <- function(x){
  loc.bt <- gev.param[1:(nk)]
  scale.bt <- exp(gev.param[(nk+1):(2*nk)])
  shape.bt <- gev.param[(2*nk+1):(3*nk)]
  
  rl.fun50 <- 0
  for (i in 1:nk){
    rl.fun50 <- rl.fun50 + HMMest$delta[i] * pevd(x, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                  type = c("GEV"),npy = 12)
  }
  (rl.fun50 - 1+ p50)
}


rl.100 <- function(x){
  loc.bt <- gev.param[1:(nk)]
  scale.bt <- exp(gev.param[(nk+1):(2*nk)])
  shape.bt <- gev.param[(2*nk+1):(3*nk)]
  
  rl.fun100 <- 0
  for (i in 1:nk){
    rl.fun100 <- rl.fun100 + HMMest$delta[i] * pevd(x, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                    type = c("GEV"),npy = 12)
  }
  (rl.fun100 - 1+ p100)
}


rl.200 <- function(x){
  loc.bt <- gev.param[1:(nk)]
  scale.bt <- exp(gev.param[(nk+1):(2*nk)])
  shape.bt <- gev.param[(2*nk+1):(3*nk)]
  
  rl.fun200 <- 0
  for (i in 1:nk){
    rl.fun200 <- rl.fun200 + HMMest$delta[i] * pevd(x, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                    type = c("GEV"),npy = 12)
  }
  (rl.fun200 - 1+ p200)
}


rl.500 <- function(x){
  loc.bt <- gev.param[1:(nk)]
  scale.bt <- exp(gev.param[(nk+1):(2*nk)])
  shape.bt <- gev.param[(2*nk+1):(3*nk)]
  
  rl.fun500 <- 0
  for (i in 1:nk){
    rl.fun500 <- rl.fun500 + HMMest$delta[i] * pevd(x, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                    type = c("GEV"),npy = 12)
  }
  (rl.fun500 - 1+ p500)
}

uniroot(rl.50,c(0,5000))$root
uniroot(rl.100,c(0,8000))$root
uniroot(rl.200,c(0,10000))$root
uniroot(rl.500,c(0,15000))$root





## Bootstrap return levels
B <- 10000
p50 <- 1/(50*12)
p100 <- 1/(100*12)
p200 <- 1/(200*12)
p500 <- 1/(500*12)

ret.l.50 <- rep(0,B)
ret.l.100 <- rep(0,B)
ret.l.200 <- rep(0,B)
ret.l.500 <- rep(0,B)
for (b in 1:B){
  gev.param.bt <- rmnorm(n = 1, mean = gev.param, varcov=gev.varc) 
  loc.bt <- gev.param.bt[1:(nk)]
  scale.bt <- exp(gev.param.bt[(nk+1):(2*nk)])
  shape.bt <- gev.param.bt[(2*nk+1):(3*nk)]
  
  rl.50 <- function(x){
    loc.bt <- gev.param.bt[1:(nk)]
    scale.bt <- exp(gev.param.bt[(nk+1):(2*nk)])
    shape.bt <- gev.param.bt[(2*nk+1):(3*nk)]
    
    rl.fun50 <- 0
    for (i in 1:nk){
      rl.fun50 <- rl.fun50 + HMMest$delta[i] * pevd(x, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                    type = c("GEV"),npy = 12)
    }
    (rl.fun50 - 1+ p50)
  }
  
  ret.l.50[b] <- uniroot(rl.50,c(0,100000))$root
  
  
  rl.100 <- function(x){
    loc.bt <- gev.param.bt[1:(nk)]
    scale.bt <- exp(gev.param.bt[(nk+1):(2*nk)])
    shape.bt <- gev.param.bt[(2*nk+1):(3*nk)]
    
    rl.fun100 <- 0
    for (i in 1:nk){
      rl.fun100 <- rl.fun100 + HMMest$delta[i] * pevd(x, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                      type = c("GEV"),npy = 12)
    }
    (rl.fun100 - 1+ p100)
  }
  
  ret.l.100[b] <- uniroot(rl.100,c(0,500000))$root
  
  
  rl.200 <- function(x){
    loc.bt <- gev.param.bt[1:(nk)]
    scale.bt <- exp(gev.param.bt[(nk+1):(2*nk)])
    shape.bt <- gev.param.bt[(2*nk+1):(3*nk)]
    
    rl.fun200 <- 0
    for (i in 1:nk){
      rl.fun200 <- rl.fun200 + HMMest$delta[i] * pevd(x, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                      type = c("GEV"),npy = 12)
    }
    (rl.fun200 - 1+ p200)
  }
  
  ret.l.200[b] <- uniroot(rl.200,c(0,1000000))$root
  
  rl.500 <- function(x){
    loc.bt <- gev.param.bt[1:(nk)]
    scale.bt <- exp(gev.param.bt[(nk+1):(2*nk)])
    shape.bt <- gev.param.bt[(2*nk+1):(3*nk)]
    
    rl.fun500 <- 0
    for (i in 1:nk){
      rl.fun500 <- rl.fun500 + HMMest$delta[i] * pevd(x, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                      type = c("GEV"),npy = 12)
    }
    (rl.fun500 - 1+ p500)
  }
  
  ret.l.500[b] <- uniroot(rl.500,c(0,5000000))$root
}

quantile(ret.l.50,probs = c(0.025,0.975))
quantile(ret.l.100,probs = c(0.025,0.975))
quantile(ret.l.200,probs = c(0.025,0.975))
quantile(ret.l.500,probs = c(0.025,0.975))



mean(ret.l.50); mean(ret.l.100); mean(ret.l.200); mean(ret.l.500);


#> quantile(ret.l.50,probs = c(0.025,0.975))
#2.5%     97.5% 
#  255.0531 1339.9358 
#> quantile(ret.l.100,probs = c(0.025,0.975))
#2.5%     97.5% 
#  348.4437 2479.9957 
#> quantile(ret.l.200,probs = c(0.025,0.975))
#2.5%     97.5% 
#  477.6632 4579.7366 
#> quantile(ret.l.500,probs = c(0.025,0.975))
#2.5%     97.5% 
#  721.3773 10563.9768 




source("Viterbi.hmmgev.R")
source("HMMNLM.resid.R")




## Monthly maxima

res.mth <- HMMNLM.resid(R1hrmax,HMMest,distn="gev")


#qqhmmMth <- qq(res.mth)

#save(qqhmmDay,file="qqhmmMth.image")
get(load("../HMM-GAMgevRes/qqhmmMth.image"))


postscript(paste("HMMRmthmax",nk,"states-ACF.eps",sep=""),paper="special",
           width=4.5*2/2,height=9)
par(mfrow=c(2,1), mar=c(4.2, 3, 1, 1))

res <- res.mth
qqhmm <- qqhmmMth
ylims<-range(res,qqhmm$low,qqhmm$upp)
plot(qqhmm$med,sort(res),ylab="",ylim=ylims,xlab="Expected",
     cex.axis=1.2,cex.lab=1.2,pch=16,cex=0.5)
title(ylab="Observed", line=2, cex.lab=1.2)
lines(qqhmm$med,qqhmm$med,lwd=0.6)
lines(qqhmm$med,qqhmm$low,col="red",lty=2,lwd=0.6)
lines(qqhmm$med,qqhmm$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(c) Monthly",bty="n",cex=1.5)

p.out<-length(which((sort(res)>qqhmm$upp)|(sort(res)<qqhmm$low)))/length(res)
if(p.out>0){
  message(paste("observations outside envelope (",round(100*p.out,1),"%)",sep=""))
}
if(p.out==0){ message("no observations outside envelope") }

inc<-1.5
lag.max <- 50
acf.hmm <- acf(res,lag.max=lag.max,plot=F)
acf.hmm.ci <- acf.ci(res)
ylims<-range(acf.hmm.ci,acf.hmm$acf,1)
plot(acf.hmm,main="",ci=0,ylim=ylims,cex.axis=1.2,cex.lab=1.2,lwd=0.6,
     xlab="Lag (month)",ylab="")
lines(acf.hmm.ci$low,lty=2,col="blue")
lines(acf.hmm.ci$upp,lty=2,col="blue")
title(ylab="ACF", line=2, cex.lab=1.2)
legend("topleft","(f) Monthly",bty="n",cex=1.5)
dev.off()








mc <- Viterbi.hmmgev(R1hrmax, HMM)

tt <- 1:length(R1hrmax)


mean.gev <- HMM$loct

temp <- rainbow(nk)
colsN <- temp[1:nk]
colsN[1] <- gray(0.2)
colsN <- colsN[order(mean.gev,decreasing = FALSE)]

vky <- mc$y
for (i in tt){
  for (j in 1:nk){
    if (vky[i]==j){
      vky[i] <- colsN[j]
    }}}


data <- R1hrmax



B <- 10000
ret.l.50 <- rep(0,B)
ret.l.100 <- rep(0,B)
ret.l.200 <- rep(0,B)
ret.l.500 <- rep(0,B)
loc.bt <- scale.bt <- shape.bt <- matrix(NA,10000,nk)
for (b in 1:B){
  gev.param.bt <- rmnorm(n = 1, mean = gev.param, varcov=gev.varc) 
  loc.bt[b,] <- gev.param.bt[1:(nk)]
  scale.bt[b,] <- exp(gev.param.bt[(nk+1):(2*nk)])
  shape.bt[b,] <- gev.param.bt[(2*nk+1):(3*nk)]
}

loc.bt.lo <- loc.bt.up <- NULL
scale.bt.lo <- scale.bt.up <- NULL
shape.bt.lo <- shape.bt.up <- NULL
for (i in 1:nk){
  loc.bt.lo[i] <- quantile(loc.bt[,i],0.025)
  loc.bt.up[i] <- quantile(loc.bt[,i],0.975)
  shape.bt.lo[i] <- quantile(shape.bt[,i],0.025)
  shape.bt.up[i] <- quantile(shape.bt[,i],0.975)
  scale.bt.lo[i] <- quantile(scale.bt[,i],0.025)
  scale.bt.up[i] <- quantile(scale.bt[,i],0.975)
}


loct.all <- loc.bt.l <- loc.bt.u <- NULL
scale.all <- scale.bt.l <- scale.bt.u <- NULL
shape.all <- shape.bt.l <- shape.bt.u <- NULL
for (i in 1:length(R1hrmax)){
  loct.all[i] <- HMM$loct[mc$y[i]]
  loc.bt.l[i] <- loc.bt.lo[mc$y[i]]
  loc.bt.u[i] <- loc.bt.up[mc$y[i]]
  shape.all[i] <- HMM$shape[mc$y[i]]
  shape.bt.l[i] <- shape.bt.lo[mc$y[i]]
  shape.bt.u[i] <- shape.bt.up[mc$y[i]]
  scale.all[i] <- HMM$scale[mc$y[i]]
  scale.bt.l[i] <- scale.bt.lo[mc$y[i]]
  scale.bt.u[i] <- scale.bt.up[mc$y[i]]
}


is_leap_year <- function(year) {
  if ((year %% 4 == 0 && year %% 100 != 0) || year %% 400 == 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

temp <- no.days <- wksoveryr <- NULL
for (i in 1994:2019){
  temp[i-1993] <- is_leap_year(i)
  if (temp[i-1993])
    no.days[i-1993] <- 366
  else no.days[i-1993] <- 365
  
  wksoveryr[i-1993] <- 12*(i-1993)
}


pdf(paste("HMMgevRmthmax",nk,"S-LocParam.pdf",sep=""),paper="special",
    width=4*cos (35.4/180*pi)/0.612,height=6) 
par(mfrow=c(3,1),mar=c(4.5, 4.5, 1, 1)) 
plot(1:length(R1hrmax),loct.all,type="l",ylim=c(min(loc.bt.l),max(loc.bt.u)),
     xlab="year",ylab="GEV location",axes=F)
box()
axis(1,at=wksoveryr,1994:2019)
axis(2)
lines(1:length(R1hrmax),loc.bt.l,lty="dashed",col="blue")
lines(1:length(R1hrmax),loc.bt.u,lty="dashed",col="blue")
legend("topright","(a)",bty="n")

plot(1:length(R1hrmax),scale.all,type="l",ylim=c(min(scale.bt.l),max(scale.bt.u)),
     xlab="year",ylab="GEV scale",axes=F)
box()
axis(1,at=wksoveryr,1994:2019)
axis(2)
lines(1:length(R1hrmax),scale.bt.l,lty="dashed",col="blue")
lines(1:length(R1hrmax),scale.bt.u,lty="dashed",col="blue")
legend("topright","(b)",bty="n")

plot(1:length(R1hrmax),shape.all,type="l",ylim=c(min(shape.bt.l),max(shape.bt.u)),
     xlab="year",ylab="GEV shape",axes=F)
box()
axis(1,at=wksoveryr,1994:2019)
axis(2)
lines(1:length(R1hrmax),shape.bt.l,lty="dashed",col="blue")
lines(1:length(R1hrmax),shape.bt.u,lty="dashed",col="blue")
legend("topright","(c)",bty="n")
dev.off()














## Bootstrap confidence interval of return levels at each time point conditioning on 
## the Viterbi path

library(mnormt)
gev.param <- y$estimate[1:(3*nk)]
gev.varc <- solve(y$hessian[1:(3*nk),1:(3*nk)])


mc <- Viterbi.hmmgev(R1hrmax, HMM)
tt <- 1:length(R1hrmax)

mean.gev <- HMM$loct

temp <- rainbow(nk)
colsN <- temp[1:nk]
colsN[1] <- gray(0.2)
colsN <- colsN[order(mean.gev,decreasing = FALSE)]

vky <- mc$y
mcy <- vky
for (i in tt){
  for (j in 1:nk){
    if (mcy[i]==j){
      mcy[i] <- colsN[j]
    }}}

data <- R1hrmax



# Estimated return levels at each time point
## Bootstrap return levels
B <- 10000

p50 <- 1/(50*12)
p100 <- 1/(100*12)
p200 <- 1/(200*12)
p500 <- 1/(500*12)

rl.50.state <- rl.100.state <- rl.200.state <- rl.500.state <- matrix(0,3,nk)
for (jj in 1:nk){
  rl.50 <- function(x){
    loc.bt <- gev.param[1:(nk)]
    scale.bt <- exp(gev.param[(nk+1):(2*nk)])
    shape.bt <- gev.param[(2*nk+1):(3*nk)]
    
    rl.fun50 <- pevd(x, loc = loc.bt[jj], scale = scale.bt[jj], shape = shape.bt[jj], 
                     type = c("GEV"),npy = 12)
    (rl.fun50 - 1 + p50)
  }
  
  
  rl.100 <- function(x){
    loc.bt <- gev.param[1:(nk)]
    scale.bt <- exp(gev.param[(nk+1):(2*nk)])
    shape.bt <- gev.param[(2*nk+1):(3*nk)]
    
    rl.fun100 <- pevd(x, loc = loc.bt[jj], scale = scale.bt[jj], shape = shape.bt[jj], 
                      type = c("GEV"),npy = 12)
    (rl.fun100 - 1 + p100)
  }
  
  
  rl.200 <- function(x){
    loc.bt <- gev.param[1:(nk)]
    scale.bt <- exp(gev.param[(nk+1):(2*nk)])
    shape.bt <- gev.param[(2*nk+1):(3*nk)]
    
    rl.fun200 <- pevd(x, loc = loc.bt[jj], scale = scale.bt[jj], shape = shape.bt[jj], 
                      type = c("GEV"),npy = 12)
    (rl.fun200 - 1 + p200)
  }
  
  
  rl.500 <- function(x){
    loc.bt <- gev.param[1:(nk)]
    scale.bt <- exp(gev.param[(nk+1):(2*nk)])
    shape.bt <- gev.param[(2*nk+1):(3*nk)]
    
    rl.fun500 <- pevd(x, loc = loc.bt[jj], scale = scale.bt[jj], shape = shape.bt[jj], 
                      type = c("GEV"),npy = 12)
    (rl.fun500 - 1 + p500)
  }
  
  rl.50.state[1,jj] <- uniroot(rl.50,c(0,30000000))$root
  rl.100.state[1,jj] <- uniroot(rl.100,c(0,10000000))$root
  rl.200.state[1,jj] <- uniroot(rl.200,c(0,1200000000))$root
  rl.500.state[1,jj] <- uniroot(rl.500,c(0,200000000))$root
  
  
  ret.l.50 <- rep(0,B)
  ret.l.100 <- rep(0,B)
  ret.l.200 <- rep(0,B)
  ret.l.500 <- rep(0,B)
  for (b in 1:B){
    gev.param.bt <- rmnorm(n = 1, mean = gev.param, varcov=gev.varc) 
    loc.bt <- gev.param.bt[1:(nk)]
    scale.bt <- exp(gev.param.bt[(nk+1):(2*nk)])
    shape.bt <- gev.param.bt[(2*nk+1):(3*nk)]
    
    rl.bt.50 <- function(x){
      loc.bt <- gev.param.bt[1:(nk)]
      scale.bt <- exp(gev.param.bt[(nk+1):(2*nk)])
      shape.bt <- gev.param.bt[(2*nk+1):(3*nk)]
      
      rl.fun50 <- pevd(x, loc = loc.bt[jj], scale = scale.bt[jj], shape = shape.bt[jj], 
                       type = c("GEV"),npy = 12)
      (rl.fun50 - 1 + p50)
    }
    
    ret.l.50[b] <- uniroot(rl.bt.50,c(0,300000000))$root
    
    
    rl.bt.100 <- function(x){
      loc.bt <- gev.param.bt[1:(nk)]
      scale.bt <- exp(gev.param.bt[(nk+1):(2*nk)])
      shape.bt <- gev.param.bt[(2*nk+1):(3*nk)]
      
      rl.fun100 <- pevd(x, loc = loc.bt[jj], scale = scale.bt[jj], shape = shape.bt[jj], 
                        type = c("GEV"),npy = 12)
      (rl.fun100 - 1 + p100)
    }
    
    ret.l.100[b] <- uniroot(rl.bt.100,c(0,1000000000))$root
    
    
    rl.bt.200 <- function(x){
      loc.bt <- gev.param.bt[1:(nk)]
      scale.bt <- exp(gev.param.bt[(nk+1):(2*nk)])
      shape.bt <- gev.param.bt[(2*nk+1):(3*nk)]
      
      rl.fun200 <- pevd(x, loc = loc.bt[jj], scale = scale.bt[jj], shape = shape.bt[jj], 
                        type = c("GEV"),npy = 12)
      (rl.fun200 - 1 + p200)
    }
    
    ret.l.200[b] <- uniroot(rl.bt.200,c(0,200000000000))$root
    
    rl.bt.500 <- function(x){
      loc.bt <- gev.param.bt[1:(nk)]
      scale.bt <- exp(gev.param.bt[(nk+1):(2*nk)])
      shape.bt <- gev.param.bt[(2*nk+1):(3*nk)]
      
      rl.fun500 <- pevd(x, loc = loc.bt[jj], scale = scale.bt[jj], shape = shape.bt[jj], 
                        type = c("GEV"),npy = 12)
      (rl.fun500 - 1 + p500)
    }
    
    ret.l.500[b] <- uniroot(rl.bt.500,c(0,50000000000))$root
  }
  
  rl.50.state[2:3,jj] <- quantile(ret.l.50,probs = c(0.025,0.975))
  rl.100.state[2:3,jj] <- quantile(ret.l.100,probs = c(0.025,0.975))
  rl.200.state[2:3,jj] <- quantile(ret.l.200,probs = c(0.025,0.975))
  rl.500.state[2:3,jj] <- quantile(ret.l.500,probs = c(0.025,0.975))
  
}




rl.50all <- rl.100all <- rl.200all <- rl.500all <- matrix(0,3,length(R1hrmax))
for (k in 1:length(R1hrmax)){
  rl.50all[,k] <- rl.50.state[,vky[k]]
  rl.100all[,k] <- rl.100.state[,vky[k]]
  rl.200all[,k] <- rl.200.state[,vky[k]]
  rl.500all[,k] <- rl.500.state[,vky[k]]
}



inc <- 2.5

postscript(paste("HMMmth",nk,"-ret.level.eps",sep=""),paper="special",
           height=2.25*4.5,width=4*4.5) 
par(mfrow=c(2,2))

par(mar=c(4.5, 4.5, 1.5, 1)) 
plot(1:length(R1hrmax),rl.50all[1,],type="l",ylim=c(min(rl.50all),inc*max(rl.50all)),
     xlab="year",ylab="Return level (nT/min)",axes=F,log="y",cex.lab=1.5)
box()
axis(1,at=wksoveryr[seq(1,26,by=5)],(1995:2020)[seq(1,26,by=5)],cex.axis=1.5)
axis(2,cex.axis=1.5)
lines(1:length(R1hrmax),rl.50all[2,],lty="dashed",col="black")
lines(1:length(R1hrmax),rl.50all[3,],lty="dashed",col="black")
legend("topleft","(a) HMM 50-year",bty="n",cex=1.5)

plot(1:length(R1hrmax),rl.100all[1,],type="l",ylim=c(min(rl.100all),inc*max(rl.100all)),
     xlab="year",ylab="Return level (nT/min)",axes=F,log="y",cex.lab=1.5)
box()
axis(1,at=wksoveryr[seq(1,26,by=5)],(1995:2020)[seq(1,26,by=5)],cex.axis=1.5)
axis(2,cex.axis=1.5)
lines(1:length(R1hrmax),rl.100all[2,],lty="dashed",col="black")
lines(1:length(R1hrmax),rl.100all[3,],lty="dashed",col="black")
legend("topleft","(b) HMM 100-year",bty="n",cex=1.5)

plot(1:length(R1hrmax),rl.200all[1,],type="l",ylim=c(min(rl.200all),inc*max(rl.200all)),
     xlab="year",ylab="Return level (nT/min)",axes=F,log="y",cex.lab=1.5)
box()
axis(1,at=wksoveryr[seq(1,26,by=5)],(1995:2020)[seq(1,26,by=5)],cex.axis=1.5)
axis(2,cex.axis=1.5)
lines(1:length(R1hrmax),rl.200all[2,],lty="dashed",col="black")
lines(1:length(R1hrmax),rl.200all[3,],lty="dashed",col="black")
legend("topleft","(c) HMM 200-year",bty="n",cex=1.5)

plot(1:length(R1hrmax),rl.500all[1,],type="l",ylim=c(min(rl.500all),inc*max(rl.500all)),
     xlab="year",ylab="Return level (nT/min)",axes=F,log="y",cex.lab=1.5)
box()
axis(1,at=wksoveryr[seq(1,26,by=5)],(1995:2020)[seq(1,26,by=5)],cex.axis=1.5)
axis(2,cex.axis=1.5)
#axis(2,at=log(c(20000)),20000,cex.axis=1.5)
lines(1:length(R1hrmax),rl.500all[2,],lty="dashed",col="black")
lines(1:length(R1hrmax),rl.500all[3,],lty="dashed",col="black")
legend("topleft","(d) HMM 500-year",bty="n",cex=1.5)
dev.off()














## Supplement plots for GAM
load("../HMM-GAMgevRes/results.RData")

# diagnostics for the rest


postscript(paste("figureS2-gamQQ-ACF.eps",sep=""),paper="special",
           width=4.5*2/2*3,height=9)
par(mfcol=c(2,3), mar=c(4.2, 3, 1, 1))

ylims<-range(pr1,qq1$low,qq1$upp)

plot(qq1$med,sort(pr1),ylab="",ylim=ylims,xlab="Expected",cex.axis=1.2,cex.lab=1.2,
     pch=16,cex=0.5)
title(ylab="Observed", line=2, cex.lab=1.2)
lines(qq1$med,qq1$med,lwd=0.6)
lines(qq1$med,qq1$low,col="red",lty=2,lwd=0.6)
lines(qq1$med,qq1$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(a) Hourly",bty="n",cex=1.5)

p.out<-length(which((sort(pr1)>qq1$upp)|(sort(pr1)<qq1$low)))/n1
if(p.out>0){
  message(paste("observations outside envelope (",round(100*p.out,1),"%)",sep=""))
}
if(p.out==0){ message("no observations outside envelope") }

inc<-1.5
lag.max<-50
acf1<-acf(pr1,lag.max=lag.max,plot=F)
ylims<-range(acf.ci1,1)
plot(acf1,main="",ci=0,ylim=ylims,cex.axis=1.2,cex.lab=1.2,lwd=0.6,
     xlab="Lag (hour)",ylab="")
lines(acf.ci1$low,lty=2,col="blue")
lines(acf.ci1$upp,lty=2,col="blue")
title(ylab="ACF", line=2, cex.lab=1.2)
legend("topleft","(d) Hourly",bty="n",cex=1.5)

ylims<-range(pr2,qq2$low,qq2$upp)

plot(qq2$med,sort(pr2),ylab="",ylim=ylims,xlab="Expected",cex.axis=1.2,cex.lab=1.2,
     pch=16,cex=0.5)
title(ylab="Observed", line=2, cex.lab=1.2)
lines(qq2$med,qq2$med,lwd=0.6)
lines(qq2$med,qq2$low,col="red",lty=2,lwd=0.6)
lines(qq2$med,qq2$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(b) Daily",bty="n",cex=1.5)

p.out<-length(which((sort(pr2)>qq2$upp)|(sort(pr2)<qq2$low)))/n2
if(p.out>0){
  message(paste("observations outside envelope (",round(100*p.out,1),"%)",sep=""))
}
if(p.out==0){ message("no observations outside envelope") }

inc<-1.5
lag.max<-50
acf2<-acf(pr2,lag.max=lag.max,plot=F)
ylims<-range(acf.ci2,1)
plot(acf2,main="",ci=0,ylim=ylims,cex.axis=1.2,cex.lab=1.2,lwd=0.6,
     xlab="Lag (day)",ylab="")
lines(acf.ci2$low,lty=2,col="blue")
lines(acf.ci2$upp,lty=2,col="blue")
title(ylab="ACF", line=2, cex.lab=1.2)
legend("topleft","(e) Daily",bty="n",cex=1.5)

ylims<-range(pr3,qq3$low,qq3$upp)

plot(qq3$med,sort(pr3),ylab="",ylim=ylims,xlab="Expected",cex.axis=1.2,cex.lab=1.2,
     pch=16,cex=0.5)
title(ylab="Observed", line=2, cex.lab=1.2)
lines(qq3$med,qq3$med,lwd=0.6)
lines(qq3$med,qq3$low,col="red",lty=2,lwd=0.6)
lines(qq3$med,qq3$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(c) Weekly",bty="n",cex=1.5)

p.out<-length(which((sort(pr3)>qq3$upp)|(sort(pr3)<qq3$low)))/n3
if(p.out>0){
  message(paste("observations outside envelope (",round(100*p.out,1),"%)",sep=""))
}
if(p.out==0){ message("no observations outside envelope") }

inc<-1.5
lag.max<-50
acf3<-acf(pr3,lag.max=lag.max,plot=F)
ylims<-range(acf.ci3,1)
plot(acf3,main="",ci=0,ylim=ylims,cex.axis=1.2,cex.lab=1.2,lwd=0.6,
     xlab="Lag (day)",ylab="")
lines(acf.ci3$low,lty=2,col="blue")
lines(acf.ci3$upp,lty=2,col="blue")
title(ylab="ACF", line=2, cex.lab=1.2)
legend("topleft","(f) Weekly",bty="n",cex=1.5)

dev.off()











