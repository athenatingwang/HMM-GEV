library(HiddenMarkov)
library(extRemes)
library(lubridate)


## Load maxima data
load(file="../MaximaData/R1wkmax.image")
R1hrmax <- R1wkmax
len.data <- length(R1hrmax)



HMM <- get(load("../HMM-GAMgevRes/HMMgevWk4est378.image"))

## Rerun HMMgev to get a better convergence
source("hmm.evd.R")
nk <- 4
Pit <- log(HMM$Pi/diag(HMM$Pi))
Pi0 <- as.vector(Pit[!diag(nk)])
loct <- HMM$loct-0.05
shapet <- HMM$shape+0.05
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
  gev.param.bt <- rmnorm(n = 1, mean = gev.param, varcov=gev.varc) 
  loc.bt <- gev.param.bt[1:(nk)]
  scale.bt <- exp(gev.param.bt[(nk+1):(2*nk)])
  shape.bt <- gev.param.bt[(2*nk+1):(3*nk)]
  
  for (i in 1:nk){
    ret.l.50[b] <- ret.l.50[b] + HMMest$delta[i] * rlevd(50*365.25/7, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                           type = c("GEV"),npy = 365.25/7)
    ret.l.100[b] <- ret.l.100[b] + HMMest$delta[i] * rlevd(100*365.25/7, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                             type = c("GEV"),npy = 365.25/7)
    ret.l.200[b] <- ret.l.200[b] + HMMest$delta[i] * rlevd(200*365.25/7, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                             type = c("GEV"),npy = 365.25/7)
    ret.l.500[b] <- ret.l.500[b] + HMMest$delta[i] * rlevd(500*365.25/7, loc = loc.bt[i], scale = scale.bt[i], shape = shape.bt[i], 
                                                             type = c("GEV"),npy = 365.25/7)
  }
}

quantile(ret.l.50,probs = c(0.025,0.975))
quantile(ret.l.100,probs = c(0.025,0.975))
quantile(ret.l.200,probs = c(0.025,0.975))
quantile(ret.l.500,probs = c(0.025,0.975))

mean(ret.l.50); mean(ret.l.100); mean(ret.l.200); mean(ret.l.500);


#> quantile(ret.l.50,probs = c(0.025,0.975))
#2.5%    97.5% 
#  277.4330 798.1694 
#> quantile(ret.l.100,probs = c(0.025,0.975))
#2.5%     97.5% 
#  388.0233 1277.0890 
#> quantile(ret.l.200,probs = c(0.025,0.975))
#2.5%    97.5% 
#  545.407 2044.240 
#> quantile(ret.l.500,probs = c(0.025,0.975))
#2.5%     97.5% 
#  856.4864 3860.7847 
#> 
#  > mean(ret.l.50); mean(ret.l.100); mean(ret.l.200); mean(ret.l.500);
#[1] 475.2183
#[1] 715.459
#[1] 1082.106
#[1] 1881.667



source("Viterbi.hmmgev.R")
source("HMMNLM.resid.R")


res <- HMMNLM.resid(R1hrmax,HMMest,distn="gev")



load("../HMM-GAMgevRes/GAMresults.RData")


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





postscript(paste("HMMRwkmax",nk,"states-GAMACF.eps",sep=""),paper="special",
           width=3*2,height=5)
par(mfrow=c(2,2), mar=c(4.4, 4.4, 1.5, 1))

qqhmm<-qq(res)
ylims<-range(res,qqhmm$low,qqhmm$upp)
plot(qqhmm$med,sort(res),ylab="Observed",ylim=ylims,xlab="Expected",
     cex.axis=1.2,cex.lab=1.2,pch=16,cex=0.5,col=ycolor)
lines(qqhmm$med,qqhmm$med,lwd=0.6)
lines(qqhmm$med,qqhmm$low,col="red",lty=2,lwd=0.6)
lines(qqhmm$med,qqhmm$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(a) HMM",bty="n",cex=1.5)


acf(res,lag.max=50,ci.col="red",cex.axis=1.2,cex.lab=1.2,xlab="Lag (week)")
legend("topleft","(b) HMM",bty="n",cex=1.5)
#pacf(res,ci.type="ma",cex.axis=1,cex.lab=1)
#legend("topleft","(c)",bty="n")


ylims<-range(pr4,qq4$low,qq4$upp)

plot(qq4$med,sort(pr4),ylab="Observed",ylim=ylims,xlab="Expected",
     cex.axis=1.2,cex.lab=1.2,pch=16,cex=0.5)
lines(qq4$med,qq4$med,lwd=0.6)
lines(qq4$med,qq4$low,col="red",lty=2,lwd=0.6)
lines(qq4$med,qq4$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(c) GAM",bty="n",cex=1.5)

acf(pr4,lag.max=50,main="",cex.axis=1.2,cex.lab=1.2,ci.col="red",xlab="Lag (month)")
legend("topleft","(d) GAM",bty="n",cex=1.5)
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
for (i in 1:length(R1wkmax)){
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
  
  wksoveryr[i-1993] <- sum(no.days[1:(i-1993)])/7
}

  
postscript(paste("HMMwk",nk,"-GAMmth-S-LocParam.eps",sep=""),paper="special",
    width=12,height=8) 
par(mfrow=c(3,3))

par(mar=c(4.5, 4.5, 1.5, 1)) 
plot(1:length(R1wkmax),loct.all,type="l",ylim=c(min(loc.bt.l),15),
     xlab="year",ylab="GEV location",axes=F,cex.lab=1.5)
box()
axis(1,at=wksoveryr[seq(1,26,by=5)],(1995:2020)[seq(1,26,by=5)],cex.axis=1.5)
axis(2,cex.axis=1.5)
lines(1:length(R1wkmax),loc.bt.l,lty="dashed",col="black")
lines(1:length(R1wkmax),loc.bt.u,lty="dashed",col="black")
legend("topright","(a) HMM",bty="n",cex=1.5)

plot(1:length(R1wkmax),scale.all,type="l",ylim=c(min(scale.bt.l),9),
     xlab="year",ylab="GEV scale",axes=F,cex.lab=1.5)
box()
axis(1,at=wksoveryr[seq(1,26,by=5)],(1995:2020)[seq(1,26,by=5)],cex.axis=1.5)
axis(2,cex.axis=1.5)
lines(1:length(R1wkmax),scale.bt.l,lty="dashed",col="black")
lines(1:length(R1wkmax),scale.bt.u,lty="dashed",col="black")
legend("topright","(b) HMM",bty="n",cex=1.5)

plot(1:length(R1wkmax),shape.all,type="l",ylim=c(min(shape.bt.l),1),
     xlab="year",ylab="GEV shape",axes=F,cex.lab=1.5)
box()
axis(1,at=wksoveryr[seq(1,26,by=5)],(1995:2020)[seq(1,26,by=5)],cex.axis=1.5)
axis(2,cex.axis=1.5)
lines(1:length(R1wkmax),shape.bt.l,lty="dashed",col="black")
lines(1:length(R1wkmax),shape.bt.u,lty="dashed",col="black")
legend("topright","(c) HMM",bty="n",cex=1.5)


load("../HMM-GAMgevRes/GAMresults.RData")

plot(m4,select=1,rug=F,scale=F,ylim=c(-1.5,1.5),ylab="GEV location",cex.axis=1.5,cex.lab=1.5)
legend("topright","(d) GAM",bty="n",cex=1.5)

plot(m4,select=3,rug=F,scale=F,trans=exp,ylim=c(0.97,1.03),ylab="GEV scale",cex.axis=1.5,cex.lab=1.5)
legend("topright","(e) GAM",bty="n",cex=1.5)

plot(m4,select=5,rug=F,scale=F,ylim=c(-0.3,0.3),ylab="GEV shape",cex.axis=1.5,cex.lab=1.5)
legend("topright","(f) GAM",bty="n",cex=1.5)

plot(m4,select=2,rug=F,scale=F,ylim=c(-15,15),ylab="GEV location",cex.axis=1.5,cex.lab=1.5)
legend("topright","(g) GAM",bty="n",cex=1.5)

plot(m4,select=4,rug=F,scale=F,trans=exp,ylim=c(0,3),ylab="GEV scale",cex.axis=1.5,cex.lab=1.5)
legend("topright","(h) GAM",bty="n",cex=1.5)

plot(m4,select=6,rug=F,scale=F,ylim=c(-0.015,0.015),ylab="GEV shape",cex.axis=1.5,cex.lab=1.5)
legend("topright","(i) GAM",bty="n",cex=1.5)

dev.off()







