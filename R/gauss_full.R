library(bestNormalize)
library(VGAM)
library(fourier)


####start####

library(mypackage)
library(bestNormalize)

sample_yy<-function(nn,x){
  nburn=0
  yy<-sqrt(1-x^2)*rnorm(1)
  for(j in 2:(nn+nburn)){
    yy[j]<-x*yy[j-1]+sqrt(1-x^2)*rnorm(1,0,1)
  }
  return(yy[(nburn+1):(nn+nburn)])
}


nn<-50
set.seed(1)
theta=.2
y <- arima.sim(model = list(order = c(1, 0, 0), ar = theta), n = nn)
lambda<-yeojohnson(y)$lambda
y_yj<-yeo.johnson(y,lambda=lambda)

hist(y)
hist(y_yj)

####start####

set.seed(1)

run_gauss<-function(iii){
theta<-runif(1)

nn<-50
nnn<-10
nburn<-50

checkR=-1
while(checkR<0){
checkR=0
checkRR<-numeric(nn)
  
  for(ii in 1:nnn){

xx<-do.call(cbind,lapply(1:1e4,function(x) as.numeric(arima.sim(model = list(order = c(1, 0, 0), ar = theta), n = nn))))
xx<-apply(xx,2,function(x) yeo.johnson(x,lambda = lambda))
checkR<-checkR+median(fr_Rm_cpp(y_yj,t(xx),seq(25,25,length.out=1)))/nnn
print("1")
}
}
(val<-log(median(checkR)))

cnt=0
theta_post<-numeric(iii)
loglik_post<-numeric(iii)
for(i in 1:iii){
  
  checkRR<-numeric(nn)
  checkR<-0
for(ii in 1:nnn){
    xx<-do.call(cbind,lapply(1:1e4,function(x) as.numeric(arima.sim(model = list(order = c(1, 0, 0), ar = theta), n = nn))))
    xx<-apply(xx,2,function(x) yeo.johnson(x,lambda = lambda))
    checkR<-checkR+median(fr_Rm_cpp(y_yj,t(xx),seq(25,25,length.out=1)))/nnn
}
  

  if(!is.na(median(checkR))&median(checkR)>0) {
  val<-log(median(checkR))
    }


  theta1<-rnorm(1,theta,sd=.1)
  if(all(theta1<1,theta1>0)){
    checkRR<-numeric(nn)
    checkR<-0
    for(ii in 1:nnn){
    xx<-do.call(cbind,lapply(1:1e4,function(x) as.numeric(arima.sim(model = list(order = c(1, 0, 0), ar = theta1), n = nn))))
    xx<-apply(xx,2,function(x) yeo.johnson(x,lambda = lambda))
      checkR<-checkR+median(fr_Rm_cpp(y_yj,t(xx),seq(25,25,length.out=1)))/nnn 

    }
      if(checkR>0) {
      val1=log(checkR)

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

return(theta_post)
}


library(parallel)
no_cores <- detectCores()-1
clust <- makeCluster(no_cores)
clusterExport(clust,c("yeojohnson","yeo.johnson","lambda",
                      "y_yj","y",
                      "fr_Rm_cpp"))
system.time(gauss_vals<-parSapply(clust, rep(300,10), run_gauss))  # 3.470099 hours
stopCluster(clust)

theta_post<-c(gauss_vals[101:300,])

par(mfrow=c(1,1))
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
hist(theta_post[501:2000],freq=F,col=c2,xlab="theta",main="summary",breaks=20,ylim=c(0,3))
abline(v=mean(theta_post[501:2000]),col="red",lwd=2)
mean(theta_post[501:(i-1)])
plot.ts(theta_post[101:(i-1)]);abline(h=mean(theta_post[101:(i-1)]),col="red",lwd=2)
plot.ts(loglik_post[101:(i-1)])

n<-50
loglikar1<-function(x){
  val1<--n/2*log(2*pi)-1/2*log((1)/(1-x^2))-(1-x^2)/(2*(1))*(y[1])^2-(n-1)/2*log((1))
  val2<-0
  for(i in 2:n) val2<-val2+(y[i]-x*y[i-1])^2
  loglik<-val1-1/(2*(1))*val2
  return(loglik)
}


theta=.2
theta_post_true<-numeric(15000)
cnt=0
for(i in 1:15000){
  theta1<-rnorm(1,theta,sd=.1)
  if(all(theta1<1,theta1>0)){
   RR<-loglikar1(theta1)-loglikar1(theta) 
  if(log(runif(1))<RR) {
    theta<-theta1
    cnt=cnt+1
  }
   }
  theta_post_true[i]<-theta
  print(paste0("iter=",i," | theta=",round(theta,3)," |accept=",round(cnt/i,3)))
}

mean(theta_post_true[2001:15000])
hist(theta_post_true[2001:15000],col=c1,add=TRUE,freq=FALSE)
abline(v=mean(theta_post_true[2001:15000]),col="blue",lwd=2)
