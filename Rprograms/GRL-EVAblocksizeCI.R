library(extRemes)
R1daymax <- get(load("../MaximaData/R1daymax.image"))
R1wkmax <- get(load("../MaximaData/R1wkmax.image"))
load(file="../MaximaData/R1.mthmax.RData")
R1mthmax <- d$rmax
R1yrmax <- get(load("../MaximaData/R1yrmax.image"))

### Daily
a <- fevd(R1daymax, method ="MLE")
a.ret <- return.level(a,return.period = c(50*365.25,100*365.25,200*365.25,500*365.25),do.ci=TRUE)
a.ret[4,]


### Weekly
b <- fevd(R1wkmax, method ="MLE")
b.ret <- return.level(b,return.period = c(50*52.14,100*52.14,200*52.14,500*52.14),do.ci=TRUE)
b.ret[4,]


### Monthly
d <- fevd(R1mthmax, method ="MLE")
d.ret <- return.level(d,return.period = c(50*12,100*12,200*12,500*12),do.ci=TRUE)
d.ret[4,]


### Yearly
e <- fevd(R1yrmax, method ="MLE")
e.ret <- return.level(e,return.period = c(50,100,200,500),do.ci=TRUE)
e.ret[4,]



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






col.fore <- c("gray","black","deepskyblue3","chocolate3")

postscript("Returnlev.blocksize.eps",paper="special",width=9*1.15,height=5*1.15,
           onefile = TRUE, horizontal = FALSE)

par(mar=c(4.5,4.5,1.5,1))

plot(1, type = "n", xlab = "Estimated return level (nT/min)",ylab = "Block size", xlim = c(-2000, 6300),
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

xleg <- rep(2900,5)
yleg <- seq(36,32,length.out=5)+1
leg <- c("50-year return level", "100-year return level", "200-year return level", "500-year return level","Median")
pchleg = c(16,16,16,16,16)
colleg = c(col.fore,"black")
for (i in 1:4)
  legend(xleg[i],yleg[i],leg[i],pch=pchleg[i],lty=1,col=colleg[i],bty="n",cex=1.5)

xleg <- 2950
yleg <- 33
leg <- c("    Median")
legend(xleg,yleg,leg,pch=16,col=1,bty="n",cex=1.5)
arrows(x0=3050, y0=30.6, 
       x1=3350, y1=30.6, code=3, angle=90, 
       length=0.05, col=1, lwd=2)
text(4750,30.6,"95% credible interval",cex=1.5)
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



