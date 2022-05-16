
set.seed(1)

#Fourier estimator 
fr_Rm<-function(x,y,R){ 
  tot=0
  for (i in 1:dim(y)[1]){
    val=NA
    for (j in 1:dim(y)[2]){
      val[j]<-sin(R[j]*(x[j]-y[i,j]))/(pi*(x[j]-y[i,j]))
    }
    tot=tot+ prod(val)
  }
  return(tot/dim(y)[1]) 
}

#or use C++ version, faster
library(fourier)


L=1e-3
yy_test<-runif(1*1e7,-L,L)
hist(yy_test,prob=TRUE)

xx_test<-0

(fr_Rm_cpp(xx_test,matrix(yy_test),20))
(true_val<-dunif(xx_test,-L,L,log = FALSE))

#for transformation
yy_alt<-matrix(yy_test)*1e4
xx_alt<-0

#checks
checkR<-errorR<-numeric(21)
rr<-seq(1,21,length.out = 21)
for(i in 1:21) {
  checkR[i]<-(fr_Rm_cpp(xx_alt,yy_alt,rr[i]))*1e4
  errorR[i]<-abs(fr_Rm_cpp(xx_alt,yy_alt,rr[i])*1e4-true_val)
  print(i/21)
}
plot.ts(errorR,xlab="R",ylab="Absolute error",xaxt='n')
axis(1, c(1, seq(1,21,length.out=21)),c(1, round(seq(1,21,length.out =21),3)),cex.axis = .6)
plot.ts(checkR,xlab="R",ylab="Density",xaxt='n')
axis(1, c(1, seq(1,21,length.out=21)),c(1, round(seq(1,21,length.out =21),3)),cex.axis = .6)
(diff_val<-median(order(abs(diff(checkR)))[1:3]))

#estimate
(fr_Rm_cpp(xx_alt,yy_alt,rr[diff_val]))*1e4

