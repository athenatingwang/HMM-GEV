### produces the envelope in the qq-plot
### as used by the code in figures.R
library(evd)

qq<-function(res){   # res = residuals
  
  pord<-function(q,p,n,j){
    p-porder(q,distn="norm",mlen=n,mean=0,sd=1,j=j,largest=F)
  }
  
  n<-length(res)
  
  med<-vector()
  low<-vector()
  upp<-vector()
  
  mult<-10
  
  int.med<-10*c(-1,1)
  int.low<-10*c(-1,0)
  int.upp<-10*c(0,1)
  
  for(j in 1:n){
    
    print(j)
    
    check<-0
    while(check==0){
      tryCatch(        {
          med[j]<-uniroot(pord,int.med,p=0.5,n=n,j=j)$root
          check<-1
        },
        error = function(e) { int.med<-int.med*mult }
        )
    }
    
    check<-0
    while(check==0){
      tryCatch(        {
        low[j]<-uniroot(pord,int.med,p=0.025,n=n,j=j)$root
        check<-1
      },
      error = function(e) { int.med<-int.med*mult }
      )
    }
    
    check<-0
    while(check==0){
      tryCatch(        {
        upp[j]<-uniroot(pord,int.med,p=0.975,n=n,j=j)$root
        check<-1
      },
      error = function(e) { int.med<-int.med*mult }
      )
    }
    
  }
  
  # output is the 50th, 2.5th and 97.5th percentiles
  # with the 2.5th and 97.5th providing the envelope
  
  list(med=med,low=low,upp=upp)
  
}
