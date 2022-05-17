library(matrixStats)
library(bestNormalize)
library(VGAM)
library(fourier)

set.seed(1)
ii<-50
a=1
b=1
p=.1
y<-rbinom(ii,5,p)

y<-y+rnorm(ii,0,.1)

lambda<-yeojohnson(y)$lambda
y_yj<-yeo.johnson(y,lambda)
hist(y_yj)

curve(dbeta(x,a+sum(y),b+5*ii-sum(y)))
abline(v=(a+sum(y))/(a+sum(y)+b+5*ii-sum(y)),col="blue")


####start####

set.seed(1)
theta<-runif(1)
nn<-10

val<--Inf
while(val==-Inf||is.na(val)){
  checkR<-numeric(ii)
  for(iii in 1:nn){
xx<-matrix(rbinom(ii*1e4,5,theta),ncol = ii)
xx<-xx+rnorm(length(c(xx)),0,0.1)
xx<-(apply(xx,1,yeo.johnson,lambda=lambda))

# #if indep
for(j in 1:ii) {
  checkR_temp<-fr_R_cpp(y_yj[j],as.matrix(t(xx)[,j]),seq(25,25,length.out=1))
  checkR[j]<-checkR[j]+checkR_temp
  
}
}
(val<-sum(log((checkR)/nn)))
print("1") 
}

cnt=0
theta_post<-numeric(2000)
loglik_post<-numeric(2000)
for(i in 1:2000){
  checkR<-numeric(ii)
  for(iii in 1:nn){
    xx<-matrix(rbinom(ii*1e4,5,theta),ncol = ii)
    xx<-xx+rnorm(length(c(xx)),0,.01)
    xx<-(apply(xx,1,yeo.johnson,lambda=lambda))
    # #if indep
    for(j in 1:ii) {
      checkR_temp<-fr_R_cpp(y_yj[j],as.matrix(t(xx)[,j]),seq(25,25,length.out=1))
      checkR[j]<-checkR[j]+checkR_temp
    }

  }
  checkR<-sum(log(checkR/nn))
  if(checkR!=-Inf&!is.na(median(checkR))){
    val<-checkR
  }

  theta1<-rnorm(1,theta,sd=.1)
  if(all(theta1<1,theta1>0)){
    checkR<-numeric(ii)
    for(iii in 1:nn){
    xx<-matrix(rbinom(ii*1e4,5,theta1),ncol = ii)
    xx<-xx+rnorm(length(c(xx)),0,0.1)
    xx<-(apply(xx,1,yeo.johnson,lambda=lambda))
    # #if indep
    for(j in 1:ii) {
      checkR_temp<-fr_R_cpp(y_yj[j],as.matrix(t(xx)[,j]),seq(25,25,length.out=1))
      checkR[j]<-checkR[j]+checkR_temp
    }
    }
    
    checkR<-sum(log(checkR/nn))
    val1<- NA
    if(checkR!=-Inf&!is.na(median(checkR))){
      val1<-median(checkR)
    }
    if(!is.na(val1)){

    RR<-(val1-val) 
    
    if(log(runif(1))<RR) {
      val=val1
      theta=theta1
      cnt=cnt+1
      
      
    }
    
  }else{
    val1=val
    theta1=theta
  }
  }
  theta_post[i]<-theta
  loglik_post[i]<-val  
  
  print(paste0("iter=",i," | theta=",round(theta,3)," |accept=",round(cnt/i,3)))
}

par(mfrow=c(1,1))
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hist(theta_post[1:(i-1)],freq=F,col=c2,xlab=expression(theta),main="",ylim=c(0,20))
curve(dbeta(x,a+sum(y),b+5*ii-sum(y)),add=TRUE,lwd=2,col="blue")
abline(v=mean(theta_post[1:(i-1)]),col="red",lwd=2)
abline(v=(a+sum(y))/(a+sum(y)+b+5*ii-sum(y)),lwd=2,col="blue")

plot.ts(theta_post[1:(i-1)])
plot.ts(loglik_post[1:(i-1)])
