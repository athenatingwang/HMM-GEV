### produces the envelope in the acf-plot
### as used by the code in figures.R

acf.ci<-function(res,lag.max=50){   # res = residuals
  
  n0<-length(res)
  
  n<-n0-(0:lag.max)
  df<-n-2
  
  t.low<-qt(0.025,df)
  t.upp<-qt(0.975,df)
  
  low<-t.low/sqrt(df+t.low^2)
  upp<-t.upp/sqrt(df+t.upp^2)
  
  # output is the lower and upper line of the envelope
  
  list(low=low,upp=upp)
  
}