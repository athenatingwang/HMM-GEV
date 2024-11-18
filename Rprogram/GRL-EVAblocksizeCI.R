setwd("B:/Research/MBIE2020CraigRodger/RprogramEVA")

data <- eyrx <- eyry <- R1 <- list()
R1max <- NULL
eyrx.all <- eyry.all <- R1.all <- NULL

count <- 1
for(i in 1994:2019){
  data[[count]] <- read.csv(paste0("../Data/DataMissingImputed/eyrxyz", i, "_fil.csv"))
  eyrx[[count]] <- data[[count]]$eyrx
  eyry[[count]] <- data[[count]]$eyry

  eyrx.all <- append(eyrx.all,eyrx[[count]])
  eyry.all <- append(eyry.all,eyry[[count]])
  
  count <- count+1
}
R1.all <- sqrt(diff(eyrx.all)^2+diff(eyry.all)^2)

R1daymax <- NULL
for (i in 1:floor(length(R1.all)/60/24)){
  R1daymax[i] <- max(R1.all[((i-1)*60*24+1):(i*60*24)])
}

R1wkmax <- NULL
for (i in 1:floor(length(R1.all)/60/24/7)){
  R1wkmax[i] <- max(R1.all[((i-1)*60*24*7+1):(i*60*24*7)])
}

R12wkmax <- NULL
for (i in 1:floor(length(R1.all)/60/24/14)){
  R12wkmax[i] <- max(R1.all[((i-1)*60*24*14+1):(i*60*24*14)])
}

R1mthmax <- NULL
for (i in 1:floor(length(R1.all)/60/24/30)){
  R1mthmax[i] <- max(R1.all[((i-1)*60*24*30+1):(i*60*24*30)])
}

R1yrmax <- NULL
for (i in 1:round(length(R1.all)/60/24/365)){
  R1yrmax[i] <- max(R1.all[((i-1)*60*24*365+1):(i*60*24*365)])
}

library(extRemes)



### Daily
a <- fevd(R1daymax, method ="MLE")
a.ret <- return.level(a,return.period = c(50*365.25,100*365.25,200*365.25,500*365.25),do.ci=TRUE)
a.ret[4,]


### Weekly
b <- fevd(R1wkmax, method ="MLE")
b.ret <- return.level(b,return.period = c(50*52.14,100*52.14,200*52.14,500*52.14),do.ci=TRUE)
b.ret[4,]


### Biweekly
c <- fevd(R12wkmax, method ="MLE")
c.ret <- return.level(c,return.period = c(50*52.14/2,100*52.14/2,200*52.14/2,500*52.14/2),do.ci=TRUE)
c.ret[4,]


### Monthly
d <- fevd(R1mthmax, method ="MLE")
d.ret <- return.level(d,return.period = c(50*12.17,100*12.17,200*12.17,500*12.17),do.ci=TRUE)
d.ret[4,]


### Yearly
e <- fevd(R1yrmax, method ="MLE")
e.ret <- return.level(e,return.period = c(50,100,200,500),do.ci=TRUE)
e.ret[4,]



qq.evd <-function(res,nsim=10^4){
  res<-sort(res)
  n<-length(res)
  res.sim<-array(NA,c(nsim,n))
  for(i in 1:nsim){ res.sim[i,]<-sort(rnorm(n)) }
  med<-apply(res.sim,2,median)
  low<-apply(res.sim,2,quantile,probs=0.025)
  upp<-apply(res.sim,2,quantile,probs=0.975)
  p.out<-length(which((res>upp)|(res<low)))/n
  if(p.out>0){
    message(paste("observations outside envelope (",round(100*p.out,1),"%)",sep=""))
  }
  else{ message("no observations outside envelope") }
  list(med=med,low=low,upp=upp)
}

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

# main code to get qq-plot for monthly maxima
# m4 = fitted model, d4 = data







col.fore <- c("gray","black","deepskyblue3","chocolate3")

postscript("Returnlev.blocksize.eps",paper="special",width=9*1.15,height=5*1.15,
           onefile = TRUE, horizontal = FALSE)

par(mar=c(4.5,4.5,1.5,1))

plot(1, type = "n", xlab = "Estimated return level (nT/min)",ylab = "Block size", xlim = c(-2000, 5300),
     ylim = c(22, 36),axes=F,cex.lab=1.5)
box()
axis(1,cex.axis=1.5)
axis(2,at=c(36,32,28,24),labels=c("day","week","month","year"),cex.axis=1.5)
#text(-2500, 36, "day",cex=1.2)

for (i in 1:4){
  points(a.ret[i,2],36-(i-1)*0.4,cex=1.2,pch=16,col=col.fore[i])
  arrows(x0=a.ret[i,1], y0=36-(i-1)*0.4, 
         x1=a.ret[i,3], y1=36-(i-1)*0.4, code=3, angle=90, 
         length=0.07, col=col.fore[i], lwd=2)
}

#text(-2500, 32, "week",cex=1.2)

for (i in 1:4){
  points(b.ret[i,2],32-(i-1)*0.4,cex=1.2,pch=16,col=col.fore[i])
  arrows(x0=b.ret[i,1], y0=32-(i-1)*0.4, 
         x1=b.ret[i,3], y1=32-(i-1)*0.4, code=3, angle=90, 
         length=0.07, col=col.fore[i], lwd=2)
}


#text(-2500, 28, "fortnight",cex=1.2)

#for (i in 1:4){
#  points(c.ret[i,2],28-(i-1)*0.4,cex=1.2,pch=16,col=col.fore[i])
#  arrows(x0=c.ret[i,1], y0=28-(i-1)*0.4, 
#         x1=c.ret[i,3], y1=28-(i-1)*0.4, code=3, angle=90, 
#         length=0.07, col=col.fore[i], lwd=2)
#}

#text(-2500, 24, "month",cex=1.2)

for (i in 1:4){
  points(d.ret[i,2],28-(i-1)*0.4,cex=1.2,pch=16,col=col.fore[i])
  arrows(x0=d.ret[i,1], y0=28-(i-1)*0.4, 
         x1=d.ret[i,3], y1=28-(i-1)*0.4, code=3, angle=90, 
         length=0.07, col=col.fore[i], lwd=2)
}

#text(-2500, 20, "year",cex=1.2)

for (i in 1:4){
  points(e.ret[i,2],24-(i-1)*0.4,cex=1.2,pch=16,col=col.fore[i])
  arrows(x0=e.ret[i,1], y0=24-(i-1)*0.4, 
         x1=e.ret[i,3], y1=24-(i-1)*0.4, code=3, angle=90, 
         length=0.07, col=col.fore[i], lwd=2)
}


#legend(x = "topright", legend = c("50-year return level", "100-year return level", "200-year return level", 
#  "500-year return level","Point estimate"),pch = c(16, 16, 16, 16, 16), col = col.fore)

xleg <- rep(2800,5)
yleg <- seq(36,32,length.out=5)+1
leg <- c("50-year return level", "100-year return level", "200-year return level", "500-year return level","Median")
pchleg = c(16,16,16,16,16)
colleg = c(col.fore,"black")
for (i in 1:4)
  legend(xleg[i],yleg[i],leg[i],pch=pchleg[i],lty=1,col=colleg[i],bty="n",cex=1.5)

xleg <- 2850
yleg <- 33
leg <- c("   Median")
legend(xleg,yleg,leg,pch=16,col=1,bty="n",cex=1.5)
arrows(x0=2900, y0=30.6, 
       x1=3200, y1=30.6, code=3, angle=90, 
       length=0.05, col=1, lwd=2)
text(4300,30.6,"95% credible interval",cex=1.5)
#legend("topleft","(a)",bty="n",cex=1.5)

dev.off()








postscript("Returnlev.blocksize.resid.eps",paper="special",width=12*0.9,
           height=10*0.9*2/4,
           onefile = TRUE, horizontal = FALSE)

par(mfrow=c(2,4))

par(mar=c(4.2,3,1,1))

pars <- a$results$par
loc <- pars["location"]
scale <- pars["scale"]
shape <- pars["shape"]

resday <- qnorm(pevd(R1daymax,loc,scale,shape,type="GEV"))
qqday<-qq(resday)
ylims<-range(resday,qqday$low,qqday$upp)
ycolor <- rep(1,length(resday))
#ycolor[qqday$ind.out] <- col.fore[3]
plot(qqday$med,sort(resday),ylab="",ylim=ylims,xlab="Expected",
     cex.axis=1.2,cex.lab=1.2,pch=16,cex=0.5,col=ycolor)
title(ylab="Observed", line=2, cex.lab=1.2)
lines(qqday$med,qqday$med,lwd=0.6)
lines(qqday$med,qqday$low,col="red",lty=2,lwd=0.6)
lines(qqday$med,qqday$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(a) daily",bty="n",cex=1.5)

acf(resday,main="",cex.axis=1.2,cex.lab=1.2,lwd=0.6,ci.col="red",xlab="Lag (day)",ylab="")
title(ylab="ACF", line=2, cex.lab=1.2)
legend("topleft","(b) daily",bty="n",cex=1.5)





pars <- b$results$par
loc <- pars["location"]
scale <- pars["scale"]
shape <- pars["shape"]

reswk <- qnorm(pevd(R1wkmax,loc,scale,shape,type="GEV"))
qqwk<-qq(reswk)
ylims<-range(reswk,qqwk$low,qqwk$upp)
ycolor <- rep(1,length(resday))
#ycolor[qqwk$ind.out] <- col.fore[3]
plot(qqwk$med,sort(reswk),ylab="",ylim=ylims,xlab="Expected",
     cex.axis=1.2,cex.lab=1.2,pch=16,cex=0.5,col=ycolor)
title(ylab="Observed", line=2, cex.lab=1.2)
lines(qqwk$med,qqwk$med,lwd=0.6)
lines(qqwk$med,qqwk$low,col="red",lty=2,lwd=0.6)
lines(qqwk$med,qqwk$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(c) weekly",bty="n",cex=1.5)

acf(reswk,main="",cex.axis=1.2,cex.lab=1.2,ci.col="red",xlab="Lag (week)",ylab="")
title(ylab="ACF", line=2, cex.lab=1.2)
legend("topleft","(d) weekly",bty="n",cex=1.5)





pars <- d$results$par
loc <- pars["location"]
scale <- pars["scale"]
shape <- pars["shape"]

resmth <- qnorm(pevd(R1mthmax,loc,scale,shape,type="GEV"))
qqmth<-qq(resmth)
ylims<-range(resmth,qqmth$low,qqmth$upp)
ycolor <- rep(1,length(resday))
#ycolor[qqmth$ind.out] <- col.fore[3]
plot(qqmth$med,sort(resmth),ylab="",ylim=ylims,xlab="Expected",
     cex.axis=1.2,cex.lab=1.2,pch=16,cex=0.5,col=ycolor)
title(ylab="Observed", line=2, cex.lab=1.2)
lines(qqmth$med,qqmth$med,lwd=0.6)
lines(qqmth$med,qqmth$low,col="red",lty=2,lwd=0.6)
lines(qqmth$med,qqmth$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(e) monthly",bty="n",cex=1.5)

acf(resmth,main="",cex.axis=1.2,cex.lab=1.2,ci.col="red",xlab="Lag (month)",ylab="")
title(ylab="ACF", line=2, cex.lab=1.2)
legend("topleft","(f) monthly",bty="n",cex=1.5)





pars <- e$results$par
loc <- pars["location"]
scale <- pars["scale"]
shape <- pars["shape"]

resyr <- qnorm(pevd(R1yrmax,loc,scale,shape,type="GEV"))
qqyr<-qq(resyr)
ylims<-range(resyr,qqyr$low,qqyr$upp)
ycolor <- rep(1,length(resday))
ycolor[qqyr$ind.out] <- col.fore[3]
plot(qqyr$med,sort(resyr),ylab="",ylim=ylims,xlab="Expected",
     cex.axis=1.2,cex.lab=1.2,pch=16,cex=0.5,col=ycolor)
title(ylab="Observed", line=2, cex.lab=1.2)
lines(qqyr$med,qqyr$med,lwd=0.6)
lines(qqyr$med,qqyr$low,col="red",lty=2,lwd=0.6)
lines(qqyr$med,qqyr$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(g) annual",bty="n",cex=1.5)

acf(resyr,main="",cex.axis=1.2,cex.lab=1.2,ci.col="red",xlab="Lag (year)",ylab="")
title(ylab="ACF", line=2, cex.lab=1.2)
legend("topleft","(h) annual",bty="n",cex=1.5)

dev.off()

































pars <- c$results$par
loc <- pars["location"]
scale <- pars["scale"]
shape <- pars["shape"]

res2wk <- qnorm(pevd(R12wkmax,loc,scale,shape,type="GEV"))
qq2wk<-qq(res2wk)
ylims<-range(res2wk,qq2wk$low,qq2wk$upp)
ycolor <- rep(1,length(resday))
#ycolor[qq2wk$ind.out] <- col.fore[3]
plot(qq2wk$med,sort(res2wk),ylab="",ylim=ylims,xlab="Expected",
     cex.axis=1.2,cex.lab=1.2,pch=16,cex=0.5,col=ycolor)
title(ylab="Observed", line=2, cex.lab=1.2)
lines(qq2wk$med,qq2wk$med,lwd=0.6)
lines(qq2wk$med,qq2wk$low,col="red",lty=2,lwd=0.6)
lines(qq2wk$med,qq2wk$upp,col="red",lty=2,lwd=0.6)
legend("topleft","(e) fortnightly",bty="n",cex=1.5)

acf(res2wk,main="",cex.axis=1.2,cex.lab=1.2,ci.col="red",xlab="Lag (fortnight)",ylab="")
title(ylab="ACF", line=2, cex.lab=1.2)
legend("topleft","(f) fortnightly",bty="n",cex=1.5)


