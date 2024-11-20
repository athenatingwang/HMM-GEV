library(HiddenMarkov)
library(extRemes)
library(lubridate)


load(file="../MaximaData/R1hrmax.image")
R1hrmax <- R1hrmax
len.data <- length(R1hrmax)



HMM <- get(load("../HMM-GAMgevRes/HMMgevHr6est780.image"))

## Rerun HMMgev to get a better convergence
source("hmm.evd.R")
nk <- 6
Pit <- log(HMM$Pi/diag(HMM$Pi))
Pi0 <- as.vector(Pit[!diag(nk)])
loct <- HMM$loct+0.02
shapet <- HMM$shape+0.02
scalet <- HMM$scale
param0<-c(loct,log(scalet),shapet,Pi0)

y <- nlm(hmm.evd,p=param0,data=R1hrmax,nk=nk,hessian = TRUE,iterlim = 10000,steptol = 1e-5,print.level=2)
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


res.hr <- HMMNLM.resid(R1hrmax,HMMest,distn="gev")
res <- res.hr


qq <- function(res,nsim=10^4){
  res<-sort(res)
  n<-length(res)
  res.sim<-array(NA,c(nsim,n))
  for(i in 1:nsim){ res.sim[i,]<-sort(rnorm(n)) }
  med<-apply(res.sim,2,median)
  low<-apply(res.sim,2,quantile,probs=0.025)
  upp<-apply(res.sim,2,quantile,probs=0.975)
  p.out<-length(which((res>upp)|(res<low)))/n
  ind.out <- which((res>upp)|(res<low))
  if(p.out>0){
    message(paste("observations outside envelope (",round(100*p.out,1),"%)",sep=""))
  }
  else{ message("no observations outside envelope") }
  list(med=med,low=low,upp=upp,ind.out=ind.out)
}




qqhmmHr <- qq(res)
qqhmm <- qqhmmHr

col.fore <- c("gray","black","deepskyblue3","chocolate3")
ycolor <- rep(1,length(res))
#ycolor[qqhmm$ind.out] <- col.fore[3]

postscript(paste("HMMRhrmax",nk,"states-ACF.eps",sep=""),paper="special",
           width=4.5*2/2,height=9)
par(mfrow=c(2,1), mar=c(4.2, 3, 1, 1))

ylims<-range(res,qqhmm$low,qqhmm$upp)
plot(qqhmm$med,sort(res),ylab="",ylim=ylims,xlab="Expected",
     cex.axis=1,cex.lab=1,pch=16,cex=0.5)
title(ylab="Observed", line=2, cex.lab=1)
lines(qqhmm$med,qqhmm$med,lwd=0.6)
lines(qqhmm$med,qqhmm$low,col="red",lty=2,lwd=0.6)
lines(qqhmm$med,qqhmm$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(a)",bty="n",cex=1.5)


acf(res,main="",lag.max=50,ci.col="red",cex.axis=1.2,cex.lab=1.2,xlab="Lag (hour)",ylab="")
title(ylab="ACF", line=2, cex.lab=1)
legend("topleft","(d)",bty="n",cex=1.5)
#pacf(res,ci.type="ma",cex.axis=1,cex.lab=1)
#legend("topleft","(c)",bty="n")
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
source("TransitionProbPlot.R")



res.day <- HMMNLM.resid(R1hrmax,HMMest,distn="gev")
res <- res.day



qqhmmDay <- qq(res)
qqhmm <- qqhmmDay

col.fore <- c("gray","black","deepskyblue3","chocolate3")
ycolor <- rep(1,length(res))
#ycolor[qqhmm$ind.out] <- col.fore[3]


postscript(paste("HMMRdaymax",nk,"states-ACF.eps",sep=""),paper="special",
           width=4.5*2/2,height=9)
par(mfrow=c(2,1), mar=c(4.2, 3, 1, 1))

ylims<-range(res,qqhmm$low,qqhmm$upp)
plot(qqhmm$med,sort(res),ylab="",ylim=ylims,xlab="Expected",
     cex.axis=1,cex.lab=1,pch=16,cex=0.5)
title(ylab="Observed", line=2, cex.lab=1)
lines(qqhmm$med,qqhmm$med,lwd=0.6)
lines(qqhmm$med,qqhmm$low,col="red",lty=2,lwd=0.6)
lines(qqhmm$med,qqhmm$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(b)",bty="n",cex=1.5)


acf(res,main="",lag.max=50,ci.col="red",cex.axis=1.2,cex.lab=1.2,xlab="Lag (day)",ylab="")
title(ylab="ACF", line=2, cex.lab=1)
legend("topleft","(e)",bty="n",cex=1.5)
#pacf(res,ci.type="ma",cex.axis=1,cex.lab=1)
#legend("topleft","(c)",bty="n")
dev.off()















####  HMM fitted to monthly maxima

## Load monthly maximum of R1 1994-2019
load(file="../MaximaData/R1.mthmax.RData")
R1hrmax <- d$rmax
len.data <- length(R1hrmax)



HMM <- get(load("../HMM-GAMgevRes/HMMgevMth3est320.image"))

## Rerun HMMgev to get a better convergence
source("hmm.evd.R")
nk <- 3
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


## Bootstrap return levels
B <- 10000
ret.l.50 <- rep(0,B)
ret.l.100 <- rep(0,B)
ret.l.200 <- rep(0,B)
ret.l.500 <- rep(0,B)
for (b in 1:B){
  gev.param.bt <- rmnorm(n = 1, mean = gev.param[1:(3*nk)], varcov=gev.varc)
  loc.bt <- gev.param.bt[1:(nk)]
  scale.bt <- exp(gev.param.bt[(nk+1):(2*nk)])
  shape.bt <- gev.param.bt[(2*nk+1):(3*nk)]
  
  for (i in 1:nk){
    ret.l.50[b] <- ret.l.50[b] + HMMest$delta[i] * rlevd(50*12, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                         type = c("GEV"),npy = 12)
    ret.l.100[b] <- ret.l.100[b] + HMMest$delta[i] * rlevd(100*12, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                           type = c("GEV"),npy = 12)
    ret.l.200[b] <- ret.l.200[b] + HMMest$delta[i] * rlevd(200*12, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                           type = c("GEV"),npy = 12)
    ret.l.500[b] <- ret.l.500[b] + HMMest$delta[i] * rlevd(500*12, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                           type = c("GEV"),npy = 12)
  }
}

quantile(ret.l.50,probs = c(0.025,0.975))
quantile(ret.l.100,probs = c(0.025,0.975))
quantile(ret.l.200,probs = c(0.025,0.975))
quantile(ret.l.500,probs = c(0.025,0.975))

mean(ret.l.50); mean(ret.l.100); mean(ret.l.200); mean(ret.l.500);


#> quantile(ret.l.50,probs = c(0.025,0.975))
#2.5%     97.5% 
#  207.0326 1226.9143 
#> quantile(ret.l.100,probs = c(0.025,0.975))
#2.5%     97.5% 
#  277.6444 2220.7086 
#> quantile(ret.l.200,probs = c(0.025,0.975))
#2.5%     97.5% 
#  370.9518 4048.7711 
#> quantile(ret.l.500,probs = c(0.025,0.975))
#2.5%     97.5% 
#  546.5566 9051.5971 
#> mean(ret.l.50); mean(ret.l.100); mean(ret.l.200); mean(ret.l.500);
#[1] 515.6376
#[1] 825.0582
#[1] 1338.624
#[1] 2590.187



source("Viterbi.hmmgev.R")
source("HMMNLM.resid.R")




## Monthly maxima

res.mth <- HMMNLM.resid(R1hrmax,HMMest,distn="gev")
res <- res.mth




qqhmmMth <- qq(res)
qqhmm <- qqhmmMth



postscript(paste("HMMRmthmax",nk,"states-ACF.eps",sep=""),paper="special",
           width=4.5*2/2,height=9)
par(mfrow=c(2,1), mar=c(4.2, 3, 1, 1))

ylims<-range(res,qqhmm$low,qqhmm$upp)
plot(qqhmm$med,sort(res),ylab="",ylim=ylims,xlab="Expected",
     cex.axis=1,cex.lab=1,pch=16,cex=0.5)
title(ylab="Observed", line=2, cex.lab=1)
lines(qqhmm$med,qqhmm$med,lwd=0.6)
lines(qqhmm$med,qqhmm$low,col="red",lty=2,lwd=0.6)
lines(qqhmm$med,qqhmm$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(c)",bty="n",cex=1.5)


acf(res,main="",lag.max=50,ci.col="red",cex.axis=1.2,cex.lab=1.2,xlab="Lag (month)",ylab="")
title(ylab="ACF", line=2, cex.lab=1)
legend("topleft","(f)",bty="n",cex=1.5)
#pacf(res,ci.type="ma",cex.axis=1,cex.lab=1)
#legend("topleft","(c)",bty="n")
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





