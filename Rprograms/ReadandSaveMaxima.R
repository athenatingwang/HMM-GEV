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


R1hrmax <- NULL
for (i in 1:floor(length(R1.all)/60)){
  R1hrmax[i] <- max(R1.all[((i-1)*60+1):(i*60)])
}
R1hrmax <- append(R1hrmax,max(R1.all[(i*60+1):(i*60+59)]))

save(R1hrmax,file="R1hrmax.image")


R1daymax <- NULL
for (i in 1:floor(length(R1.all)/60/24)){
  R1daymax[i] <- max(R1.all[((i-1)*60*24+1):(i*60*24)])
}
R1daymax <- append(R1daymax,max(R1.all[(i*60*24+1):(i*60*24+60*24-1)]))

save(R1daymax,file="R1daymax.image")


R1wkmax <- NULL
for (i in 1:floor(length(R1.all)/60/24/7)){
  R1wkmax[i] <- max(R1.all[((i-1)*60*24*7+1):(i*60*24*7)])
}
R1wkmax <- append(R1wkmax,max(R1.all[(i*60*24*7+1):(i*60*24*7+60*24*7-1)]))

save(R1wkmax,file="R1wkmax.image")


R1mthmax <- NULL
for (i in 1:floor(length(R1.all)/60/24/30)){
  R1mthmax[i] <- max(R1.all[((i-1)*60*24*30+1):(i*60*24*30)])
}
R1mthmax <- append(R1mthmax,max(R1.all[(i*60*24*30+1):(i*60*24*30+60*24*30-1)]))

save(R1mthmax,file="R1mthmax.image")


R1yrmax <- NULL
for (i in 1:round(length(R1.all)/60/24/365)){
  R1yrmax[i] <- max(R1.all[((i-1)*60*24*365+1):(i*60*24*365)])
}

save(R1yrmax,file="R1yrmax.image")

