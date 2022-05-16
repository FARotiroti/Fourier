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
par(mfrow=c(1,3))
hist(y,freq=F,breaks=50,main="",xlab="y")
y<-y+rnorm(ii,0,.1)
hist(y,freq=F,breaks=50,main="",xlab="y + noise")

lambda<-yeojohnson(y)$lambda
y_yj<-yeo.johnson(y,lambda)
hist(y_yj,freq=F,breaks=50,main="",xlab="y + noise + Yeo-Johnson")

curve(dbeta(x,a+sum(y),b+5*ii-sum(y)))
abline(v=(a+sum(y))/(a+sum(y)+b+5*ii-sum(y)),col="blue")
####start####

set.seed(1)
theta<-runif(1)
nn<-10
checkR<--1
while(checkR<0){
  
  checkR=0
  for(iii in 1:nn){
xx<-matrix(rbinom(ii*1e4,5,theta),ncol = ii)
xx<-(apply(xx,1,yeo.johnson,lambda=lambda))
checkR_temp<-median(fr_Rm_cpp(y_yj,t(xx),seq(100,100,length.out=1)))/nn
checkR<-checkR+checkR_temp
}
if(is.na(checkR)) checkR<--1

print("1") 
}
 (val<-log(median(checkR)))


cnt=0
theta_post<-numeric(5000)
loglik_post<-numeric(5000)
for(i in 1:5000){
  checkR=0
  for(iii in 1:nn){
    xx<-matrix(rbinom(ii*1e4,5,theta),ncol = ii)
    xx<-xx+rnorm(length(c(xx)),0,.1)
    xx<-(apply(xx,1,yeo.johnson,lambda=lambda))
    checkR_temp<-median(fr_Rm_cpp(y_yj,t(xx),seq(100,100,length.out=1)))/nn
    checkR<-checkR+checkR_temp
}
    if(!is.na(median(checkR))&&median(checkR)>0) {
      val<-log(median(checkR))
      }

  
  theta1<-rnorm(1,theta,sd=.1)
  if(all(theta1<1,theta1>0)){
    checkR=0
    for(iii in 1:nn){
    xx<-matrix(rbinom(ii*1e4,5,theta1),ncol = ii)
    xx<-xx+rnorm(length(c(xx)),0,0.1)
    xx<-(apply(xx,1,yeo.johnson,lambda=lambda))
    checkR_temp<-median(fr_Rm_cpp(y_yj,t(xx),seq(100,100,length.out=1)))/nn
    checkR<-checkR+checkR_temp    
    }
    
    if(is.na(median(checkR))) {
      val1<- -1
    }else{
      val1<-median(checkR)
    }
  if(val1>0) {
    val1=log(val1)
    
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




