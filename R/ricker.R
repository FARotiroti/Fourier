library(sl)
library(matrixStats)
library(VGAM)
library(bestNormalize)
library(fourier)

set.seed(3)
n.t <- 50 ## length of data run
burn.in <- 50 ## burnin length for model
n.reps <- 1e4 ## slow mixing if this lowered 

## simulate random deviates to re-use for fitting...
e <- matrix(rnorm((burn.in+n.t)*n.reps),burn.in+n.t,n.reps)
u <- runif(n.t*n.reps) 

## simulate data ....
theta <- c(3.8,0.3,10)
names(theta) <- c("r","sig.e","phi")
Y <- ricker(theta,e,u,burn.in)
y <- Y[,1]



########################################################
## The MCMC chain exploring the synthetic log likelihood
########################################################

rktrans <- NULL ## default value (stats used raw)

data(rktrans) 

set.seed(3) ## comment out for replicates!!

n.mc <- 5000 ## value for package checking
theta <- c(4,0.1,6) ## gets dynamics somewhere in right ball park
n.para <- length(theta)
th <- matrix(0,n.para,n.mc) ## storage for chain results
llr <- rep(NA,n.mc)         ## storage for log synthetic likelihood
th[,1] <- log(theta)        ## initial state
prop.sd <- c(.02,.1,.05)    ## proposal standard deviations


llr[1] <- ricker.ll(exp(th[,1]),y,e,u,burn.in=burn.in,trans=rktrans)
reject <- 0
uni <- runif(n.mc)
pname <- c("r","sig.e","phi")

for (i in 2:n.mc) { ## MCMC loop
  th[,i] <- th[,i-1] + rnorm(n.para)*prop.sd*2
  
  llr[i] <- ricker.ll(exp(th[,i]),y,e,u,burn.in,trans=rktrans)
  
  alpha <- min(1,exp(llr[i]-llr[i-1]))
  
  if (uni[i]>alpha) { ## reject
    th[,i] <- th[,i-1]
    llr[i] <- llr[i-1]
    reject <- reject + 1
  }
  
  if (i%%200==0) { ## plot the chains so far (and report acceptance rate)
    par(mfrow=c(2,2))
    start <- 1
    for (j in 1:n.para) { 
      xx <- exp(th[j,1:i])
      if (i < 5000) ylim=range(xx) else ylim=range(xx[(i/2):i])
      plot(1:i,xx,type="l",ylab=pname[j],ylim=ylim)
    }
    xx <- llr[1:i]
    if (i < 5000) ylim=range(xx) else ylim=range(xx[(i/2):i])
    plot(1:i,xx,type="l",ylab="log l_s",ylim=ylim,main=round(1-reject/i,digits=2))
  }
} ## end of the MCMC chain

1-reject/n.mc ## acceptance rate


start <- round(n.mc/2)
for (j in 1:n.para) {
  hist(exp(th[j,start:n.mc]),xlab=pname[j]);print(mean(exp(th[j,start:n.mc])))
  abline(v=mean(exp(th[j,start:n.mc])),col="blue",lwd=2)
}
hist(llr[start:n.mc],xlab="log l_s");mean(llr[start:n.mc])
abline(v=mean(llr[start:n.mc]),col="blue",lwd=2)


#fourier

#mcmc

set.seed(1)
y<-y+rnorm(length(y),0,sd=.1)
lambda<-yeojohnson(y)$lambda
y_yj<-yeo.johnson(y,lambda)

NN<-2000
samp_mcmc<-function(NN){
  nn<-10
  burn=0
  theta <-  c(3.8,.3,10)
  # theta[1]<-runif(1,3,4)
  # theta[2]<-runif(1,.1,.4)
  # theta[3]<-runif(1,7,9)
  

  
  llr <- rep(NA,NN)         ## storage for log synthetic likelihood
  th <- matrix(0,3,NN) ## storage for chain results
  se <- rep(NA,NN)  
  th[,1] <- log(theta)        ## initial state
  
      val1=-Inf
      while(val1==-Inf||is.na(median(val1))){
        checkRR<-numeric(n.t)
for(ii in 1:nn){
    e <- matrix(rnorm((burn.in+n.t)*n.reps),burn.in+n.t,n.reps)
    u <- runif(n.t*n.reps)
    Y <- ricker(exp(th[,1]), e, u, burn.in)

    Y=Y+rnorm(length(Y),0,sd = 0.1)
    Y<-apply(Y,2,function(x) yeo.johnson(x,lambda = lambda))
  
  checkRR[1]<-checkRR[1]+(median(fr_Rm_cpp(y_yj[1],as.matrix(t(Y)[,1]),seq(25,25,length.out=1))))
    for(j in 2:(n.t)){

      checkR_temp_n<-(median(fr_Rm_cpp(y_yj[c(j,j-1)],as.matrix(t(Y)[,c(j,j-1)]),seq(25,25,length.out=1))))
      checkR_temp_d<-(median(fr_Rm_cpp(y_yj[j-1],as.matrix(t(Y)[,j-1]),seq(25,25,length.out=1))))

      checkRR[j]=checkRR[j]+checkR_temp_n/checkR_temp_d

    }

    }

        val1<-sum(log(checkRR/nn))
        
    print("1")
  }
  
    (llr[1] <-(val1))

  prop.sd <- c(.02,.1,.05)  #on log scale
  reject <- 0
  uni <- runif(NN)
   
  for (i in 2:(NN+burn)) { ## MCMC loop
    

##prev

        checkRR_xx<--1
        checkRR<-numeric(n.t)
        for(ii in 1:nn){
        # e <- matrix(rnorm((burn.in+n.t)*n.reps),burn.in+n.t,n.reps)
        # u <- runif(n.t*n.reps)
        Y <- ricker(exp(th[,i-1]), e, u, burn.in)
        Y=Y+rnorm(length(Y),0,sd = 0.1)
        Y<-apply(Y,2,function(x) yeo.johnson(x,lambda = lambda))
        checkRR[1]<-checkRR[1]+(median(fr_Rm_cpp(y_yj[1],as.matrix(t(Y)[,1]),seq(25,25,length.out=1))))
        for(j in 2:(n.t)){

          checkR_temp_n<-(median(fr_Rm_cpp(y_yj[c(j,j-1)],as.matrix(t(Y)[,c(j,j-1)]),seq(25,25,length.out=1))))
          checkR_temp_d<-(median(fr_Rm_cpp(y_yj[j-1],as.matrix(t(Y)[,j-1]),seq(25,25,length.out=1))))

         if(checkR_temp_n<0) checkR_temp_n=0
          checkRR[j]=checkRR[j]+checkR_temp_n/checkR_temp_d

        }
        

        }

        checkR<-sum(log(checkRR/nn))
        if(checkR!=-Inf&!is.na(median(checkR))) {
          checkRR_xx<-(median(checkR))
      
        llr[i-1]= (checkRR_xx)

      }

 th[,i] <- th[,i-1]
        th_ind<-sample(1:3,1)
        th[th_ind,i] <- th[th_ind,i-1] + rnorm(1)*prop.sd[th_ind]*2
        
    val1=-1
    checkRR<-numeric(n.t)
    for(ii in 1:nn){
      # e <- matrix(rnorm((burn.in+n.t)*n.reps),burn.in+n.t,n.reps)
      # u <- runif(n.t*n.reps)
      Y <- ricker(exp(th[,i]), e, u, burn.in)
      Y=Y+rnorm(length(Y),0,sd = 0.1)
      Y<-apply(Y,2,function(x) yeo.johnson(x,lambda = lambda))
      checkRR[1]<-checkRR[1]+(median(fr_Rm_cpp(y_yj[1],as.matrix(t(Y)[,1]),seq(25,25,length.out=1))))
      for(j in 2:(n.t)){

        checkR_temp_n<-(median(fr_Rm_cpp(y_yj[c(j,j-1)],as.matrix(t(Y)[,c(j,j-1)]),seq(25,25,length.out=1))))
        checkR_temp_d<-(median(fr_Rm_cpp(y_yj[j-1],as.matrix(t(Y)[,j-1]),seq(25,25,length.out=1))))
        checkRR[j]=checkRR[j]+checkR_temp_n/checkR_temp_d

      }

    }
    
    checkR<-sum(log(checkRR/nn))
    if(checkR!=-Inf&!is.na(median(checkR))) {
      val1=(checkR)
        llr[i]=(val1)
      
      RR<-(llr[i]-llr[i-1])

      alpha <- min(1,exp(RR))
      
      if (uni[i]>alpha) { ## reject
        th[,i] <- th[,i-1]
        llr[i] <- llr[i-1]
        reject <- reject + 1
      }
      
    }else{
      th[,i] <- th[,i-1]
      llr[i] <- llr[i-1]
      reject <- reject + 1
    }
    
    
    cat(paste0("iter=",i),"|",paste0("theta=",round(exp(th[,i]),3) ), "|",   
        paste0("accept=",round(1-reject/i,3)),"\n")
    
    
  }
  
  #   if(j%%1e2==0) cat(paste0("iter=",j),"|",paste0("cnt=",cnt),"\n")
  etime<-Sys.time() 
  return(list(theta=exp(th),llr=llr))
}

par(mfrow=c(2,2))
1-reject/(i-1) ## acceptance rate
plot.ts(exp(th[1,1:(i-1)]));mean((exp(th[1,1:(i-1)])));abline(h=mean((exp(th[1,1:(i-1)]))),col="blue",lwd=2)
plot.ts(exp(th[2,1:(i-1)]));mean((exp(th[2,1:(i-1)])));abline(h=mean((exp(th[2,1:(i-1)]))),col="blue",lwd=2)
plot.ts(exp(th[3,1:(i-1)]));mean((exp(th[3,1:(i-1)])));abline(h=mean((exp(th[3,1:(i-1)]))),col="blue",lwd=2)

plot.ts(llr[1:(i-1)])
exp(th[,order((llr[1:(i-1)]),decreasing = TRUE)])[,1:10]
exp(th[,which.max((llr[1:(i-1)]))])


#parallel
library(parallel)
no_cores <- detectCores()-1
clust <- makeCluster(no_cores)
clusterExport(clust,c("ricker","yeojohnson","yeo.johnson","lambda",
                      "n.reps", "burn.in","fr_Rm_cpp","n.t",
                      "y_yj","y"))
system.time(save_ricker<-parSapply(clust, rep(3000,10), samp_mcmc)) 

stopCluster(clust)

par(mfrow=c(1,3))
hist(col=c2,do.call(c,lapply(save_ricker[1,], function(x) x[1,1001:3000])),main="",xlab = "theta1",freq=FALSE,breaks = 15,xlim=c(3.2,4.5))
abline(v=mean(do.call(c,lapply(save_ricker[1,], function(x) x[1,1001:3000]))),col="red",lwd=2)
abline(v=3.8,col="blue",lwd=2)
mean(do.call(c,lapply(save_ricker[1,], function(x) x[1,1001:3000])))
median(do.call(c,lapply(save_ricker[1,], function(x) x[1,1001:3000])))

hist(col=c2,do.call(c,lapply(save_ricker[1,], function(x) x[2,1001:3000])),main="",xlab = "theta2",freq=FALSE,breaks = 15,xlim=c(.05,0.6))
abline(v=mean(do.call(c,lapply(save_ricker[1,], function(x) x[2,1001:3000]))),col="red",lwd=2)
abline(v=0.3,col="blue",lwd=2)
mean(do.call(c,lapply(save_ricker[1,], function(x) x[2,1001:3000])))
median(do.call(c,lapply(save_ricker[1,], function(x) x[2,1001:3000])))

hist(col=c2,do.call(c,lapply(save_ricker[1,], function(x) x[3,1001:3000])),main="",xlab = "theta3",freq=FALSE,breaks = 5,xlim=c(6,14))
abline(v=mean(do.call(c,lapply(save_ricker[1,], function(x) x[3,1001:3000]))),col="red",lwd=2)
abline(v=10,col="blue",lwd=2)
mean(do.call(c,lapply(save_ricker[1,], function(x) x[3,1001:3000])))
median(do.call(c,lapply(save_ricker[1,], function(x) x[3,1001:3000])))

par(mfrow=c(1,1))
plot.ts(save_ricker[1,][[1]][1,],ylim=c(0,5))
for(i in 1:30)
lines(save_ricker[1,][[i]][1,])

par(mfrow=c(1,1))
plot.ts(save_ricker[1,][[1]][2,],ylim=c(0,1))
for(i in 1:30)
  lines(save_ricker[1,][[i]][2,])

par(mfrow=c(1,1))
plot.ts(save_ricker[1,][[1]][3,],ylim=c(5,15))
for(i in 1:30)
  lines(save_ricker[1,][[i]][3,])



max_mat<-matrix(0,10,ncol = 3)
max_vals<-numeric(10)
for(i in 1:10){
max_ind<-which.max(save_ricker[2,][[i]])
max_mat[i,]<-save_ricker[1,][[i]][,max_ind]
max_vals[i]<-save_ricker[2,][[i]][max_ind]
}
apply(max_mat,2,mean)
max_mat[which.max(max_vals),]

