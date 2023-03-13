library(fourier)
library(ks)
library(parallel)

### univariate mixture

sim_model_eval<-function(x,theta){
  dnorm(x,theta,1)/2+
    dnorm(x,theta,sd=1/10)/2
}

sim_model_grid<-function(theta){
  lvals<-log((sim_model_eval(seq(-6,6,by = .0001),theta)))
  delta <-.0001
  x <- sample(seq(-6,6,by =delta),1,replace=TRUE,prob=exp(lvals-max(lvals))) +
    runif(1,min=-delta/2,max=+delta/2) # smear out x
  return(x)
}


sim_model_n_grid<-function(theta,n){
  lvals<-log((sim_model_eval(seq(-6,6,by = .0001),theta)))
  delta <-.0001
  x <- sample(seq(-6,6,by =delta),n,replace=TRUE,prob=exp(lvals-max(lvals))) +
    runif(n,min=-delta/2,max=+delta/2) # smear out x
  return(x)
}

sim_model<-function(theta){
  w<-sample(0:1,1,replace=TRUE)
  w*rnorm(1,theta,1)+
    (1-w)*rnorm(1,theta,sd=1/10)
}


sim_model_n<-function(theta,n){
  w<-sample(0:1,n,replace=TRUE)
  w*rnorm(n,theta,1)+
    (1-w)*rnorm(n,theta,sd=1/10)
}


vals<-(sim_model_eval(seq(-6,6,by = .0001),0))
sum(vals)*.0001
hist(sample(seq(-6,6,by = .0001),5e5,prob = vals,replace = TRUE),freq=F,breaks = 100)
lines(seq(-6,6,.0001),sim_model_eval(seq(-6,6,.0001),0),lwd=3)

theta<-0
samp_data<-sim_model_n(theta,1e4)
hist(samp_data,freq=F)

# samp_data[samp_data==0]=1e-12
check_vals<-numeric(20)
for(j in 1:20) check_vals[j]<- fr_Rm_cpp(0,as.matrix(samp_data)*10,j)*10
plot.ts(check_vals)
abline(h=sim_model_eval(0,theta),col="blue")
abline(h=median(fr_Rm_cpp(0,as.matrix(samp_data)*10,5:20)*10))
abline(h=median(fr_Rm_cpp(0,as.matrix(samp_data)*10,5:5)*10))
for(j in 1:20) check_vals[j]<- fr_Rm_cpp((0-mean(samp_data))/sd(samp_data),as.matrix(xx)*10,j)*10/sd(samp_data)
plot.ts(check_vals)
abline(h=sim_model_eval(0,0))

(val<-sim_model_eval(0,0))
fr_Rm_cpp(0,as.matrix(samp_data)*10,10)*10


no.cores<-detectCores()-1
cl <- makeCluster(getOption("cl.cores", no.cores))

clusterExport(cl,c("sim_model_n"))

invisible(clusterEvalQ(cl, library(fourier)))
invisible(clusterEvalQ(cl, library(VGAM)))

f_est<-function(theta_1){
 xx<- sim_model_n(theta_1,1e5)
return(median(fr_Rm_cpp(0,as.matrix(xx)*10,5:5))*10)
}

val<-mean(unlist(parSapply(cl,rep(theta,10),f_est)))


xx<-runif(1e6,-10,10)
xx_d<-abs(xx-0)
xx_q<-quantile(xx_d,0.01/100)
xx_q<-0.2


cnt=0

post_samp_f<-0
# post_samp_abc<-numeric(1e5)
# post_samp_kde<-numeric(1e5)

for(j in 1:1e5){
  theta_1<-rnorm(1,theta,.5)
  
  # val_1=sim_model_eval(0,theta_1)
  
  xx<-sim_model_n(theta_1,1e4)

  val_1<- mean(unlist(parSapply(cl,rep(theta_1,10),f_est)))
  # val_1<- kde(x=xx,eval.points = 0)$estimate
  
if(runif(1)<val_1/val){
# if( abs(mean(xx)-0)<xx_q){
val=val_1
theta=theta_1
cnt=cnt+1
  }
  post_samp_f[j]<-theta
  if(j%%100==0) cat(paste0("iter=",j),"|",paste0("cnt=",cnt),"\n")
}

stopCluster(cl)


### multivariate mixture

library(mvtnorm)

sim_model_eval<-function(x,theta){
  dmvnorm(x,theta,diag(length(theta)))/2+
    dmvnorm(x,theta,0.01*diag(length(theta)))/2
}

xy<-seq(-3,3,by = .01)
z<-matrix(0,length(xy),length(xy))
for(i in 1:length(xy)){
  for(j in 1:length(xy)){
    z[i,j]<-sim_model_eval(rep(0,p),c(xy[i],xy[j]))
    
  }
}


sim_model<-function(theta){
  w<-sample(0:1,1,replace=TRUE)
  w*rmvnorm(1,theta,diag(length(theta)))+
    (1-w)*rmvnorm(1,theta,0.01*diag(length(theta)))
}


sim_model_n<-function(theta,n){
  w<-sample(0:1,n,replace=TRUE)
  w*rmvnorm(n,theta,diag(length(theta)))+
    (1-w)*rmvnorm(n,theta,0.01*diag(length(theta)))
}


p<-2
theta<-rep(0,p)
samp_data<-sim_model_n(theta,1e5)
# plot(samp_data)

check_vals<-numeric(20)
for(j in 1:20) check_vals[j]<- fr_Rm_cpp(rep(0,p),as.matrix(samp_data)*10,j)*10^p
plot.ts(check_vals)
abline(h=sim_model_eval(rep(0,2),theta),col="blue")

(val<-sim_model_eval(rep(0,p),theta))
median(fr_Rm_cpp(rep(0,p),as.matrix(samp_data)*10,5:15)*10^p)
kde(x=samp_data,eval.points = rep(0,p))$estimate



no.cores<-detectCores()-1
cl <- makeCluster(getOption("cl.cores", no.cores))

clusterExport(cl,c("sim_model_n"))

invisible(clusterEvalQ(cl, library(fourier)))
invisible(clusterEvalQ(cl, library(VGAM)))
invisible(clusterEvalQ(cl, library(mvtnorm)))
invisible(clusterEvalQ(cl, library(ks)))

f_est<-function(theta_1){
  p<-length(theta_1)
  xx<- sim_model_n(theta_1,1e4)
  return(median(fr_Rm_cpp(rep(0,p),as.matrix(xx)*10,5))*10^p)
}

f_est_kde<-function(theta_1){
  p<-length(theta_1)
  xx<- sim_model_n(theta_1,1e4)
  return(kde(x=xx,eval.points = rep(0,p))$estimate)
}

theta_list<-list()
for(i in 1:10) theta_list[[i]]<-theta
(val<-mean(unlist(parLapply(cl,theta_list,f_est))))
(val<-mean(unlist(parLapply(cl,theta_list,f_est_kde))))
sim_model_eval(rep(0,p),theta)

xx<-runif(1e6,-10,10)
xx_d<-abs(xx-0)
xx_q<-quantile(xx_d,0.01/100)
xx_q<-0.2


cnt=0
# post_samp<-matrix(0,1e4,2)
# post_samp_f<-matrix(0,1e4,2)
# post_samp_abc<-matrix(0,1e4,2)
# post_samp_kde<-matrix(0,1e4,2)

theta<-rep(0,p)

for(j in 1:1e4){
  theta_1<-rmvnorm(1,theta,.5^2*diag(p))
  
  # val_1=sim_model_eval(rep(0,p),theta_1)
  
  for(i in 1:10) theta_list[[i]]<-theta_1
  val_1<-mean(unlist(parLapply(cl,theta_list,f_est_kde)))
  # val_1<- mean(unlist(parSapply(cl,theta_list,f_est)))
  

  
  if(runif(1)<val_1/val){
    # if( abs(mean(xx)-0)<xx_q){
    val=val_1
    theta=theta_1
    cnt=cnt+1
  }
  post_samp_kde[j,]<-theta
  if(j%%100==0) cat(paste0("iter=",j),"|",paste0("cnt=",cnt),"\n")
}

stopCluster(cl)


### g-and-k

library(gk)
library(mvtnorm)
library(bestNormalize)
library(VGAM)
library(fourier)
library(parallel)


  y0<-data1000[,iii]
  hist(y0)
  
  n<-1000

muy<-mean(y0)
sdy<-sd(y0)
y<-(y0-muy)/sdy
y<-y0
hist(y)

exact=FALSE
for(jjj in 1:2){
  if(jjj==2) exact=TRUE

q_y<-quantile(y,seq(0,1,by=0.1))
q_y<-quantile(y,c(0,0.1,0.25,0.75,0.9,1))
q_y<-quantile(y,c(0,1))
q_y<-quantile(y,c(0,0.25,0.75,1))

yq_y<-NULL
for(i in 1:(length(q_y)-1)) yq_y[[i]]<-y[y>=q_y[i]&y<q_y[i+1]]
yq_y[[i]]<-c(yq_y[[i]],max(y))


lambda_y<-NULL
for(i in 1:length(yq_y)) lambda_y[[i]]<-yeojohnson(yq_y[[i]])$lambda

y_yj1_y<-
y_yj1_mu<-
y_yj1_sd<-NULL
for(i in 1:length(yq_y)){
  y_yj1_y[[i]]<-yeo.johnson(yq_y[[i]],lambda_y[[i]])
  y_yj1_mu[[i]]<-mean(y_yj1_y[[i]])
  y_yj1_sd[[i]]<-sd(y_yj1_y[[i]])
  y_yj1_y[[i]]<-(y_yj1_y[[i]]-y_yj1_mu[[i]])/y_yj1_sd[[i]]
} 

for(i in 1:length(y_yj1_y)) hist(y_yj1_y[[i]])

y_ind<-list()
y_len<-sapply(y_yj1_y,length)
cnt=1
for(i in 1:length(y_len)){  
  for(j in 1:(y_len[i])){
    y_ind[[cnt]]<-c(j,i)
    cnt=cnt+1
  }}

flambda<-function(y,lambda){
  if(y>=0) val<-abs((y+1)^(lambda-1))
  if(y<0) val<-abs((1-y)^(1-lambda))
  return(val)
}

R_y<-NULL
for(i in 1:length(yq_y)) R_y[[i]]<-floor(1+apply(as.matrix(yq_y[[i]]),1,
                                                 function(x) 10/
                                                   (y_yj1_sd[[i]]*flambda(x,lambda_y[[i]])) 
                                                 ))



f_est<-function(arg1){
  
  i<-arg1[1]
  j<-arg1[2]
  yy_yj1_y<-yeo.johnson(yy,lambda_y[[j]])
  yy_yj1_y<-(yy_yj1_y-y_yj1_mu[[j]])/y_yj1_sd[[j]]
  
  f_est1<-fr_Rm_cpp(y_yj1_y[[j]][[i]],as.matrix(yy_yj1_y),
                    seq(5,15,length.out=10))
  
  return(median(f_est1))
  
}

f_est_orig<-function(arg1){
  
  i<-arg1[1]
  
  f_est1<-fr_Rm_cpp(y[i],as.matrix(yy),5:15)
  
  return(median(f_est1))
  
}



no.cores<-detectCores()-1
# no.cores<-10
cl <- makeCluster(getOption("cl.cores", no.cores))

theta<-c(3,1,2,0.5)
a=theta[1]
b=theta[2]
g=theta[3]
k=theta[4]
c=0.8
yy<-rgk(1e4, A=a, B=b, g=g, k=k,c=c)
# yy<-yeo.johnson(yy,lambda)
# yy<-(yy-muy)/sdy

clusterExport(cl,c("y","yq_y",
                   "y_yj1_y","yy",
                   "lambda_y","R_y",
                   "flambda","muy","sdy",
                   "y_yj1_mu","y_yj1_sd"
),
  envir = environment()
)

invisible(clusterEvalQ(cl, library(fourier)))
invisible(clusterEvalQ(cl, library(VGAM)))
invisible(clusterEvalQ(cl, library(gk)))

f_est1_1<-unlist(parLapply(cl,y_ind,f_est))
# f_est1_1<-unlist(parLapply(cl,1:length(y),f_est_orig))
(f_est1_0<-which(f_est1_1<0))
f_est1_1[f_est1_0]=min(f_est1_1[f_est1_1>0])
(lval<-sum(log(f_est1_1)))

# theta<-c(3,1,2,0.5)
# a=theta[1]
# b=theta[2]
# g=theta[3]
# k=theta[4]
# sum(dgk(y0, A=a, B=b, g=g, k=k,c=c,log = TRUE))
# lval<-sum(dgk(y0, A=a, B=b, g=g, k=k,c=c,log = TRUE))


# f_est1_1<-kde(x=yy,eval.points = y)$estimate
# (f_est1_0<-which(f_est1_1<0))
# f_est1_1[f_est1_0]<-min(f_est1_1[f_est1_1>0])
# lval<-sum(log(f_est1_1))


post_samp<-matrix(0,5e4,4)
cnt=0

prior_mat<-matrix(0,4,2)
prior_mat[1,]<-c(0,10)
prior_mat[2,]<-c(0,10)
prior_mat[3,]<-c(0,10)
prior_mat[4,]<-c(0,10)

# cov_rw = matrix(c(0.00395704157562550,0.00575469238181253,-0.0428438883489085,-0.0115506423087838,
#           0.00575469238181253,0.0643904092877776,-0.000659547287101048,-0.162780899868077,
#           -0.0428438883489085,-0.000659547287101048,1.11795064638688,0.0584745415184396,
#           -0.0115506423087838,-0.162780899868077,0.0584745415184396,0.809135776947481),4,4,byrow = TRUE)

cov_rw = matrix(c(0.000324168781523849,0.00118021913416727,-0.000595390778348646,-0.00120801661101190,
                  0.00118021913416727,0.00837993122285485,0.00119368970687378,-0.00868718303162277,
                  -0.000595390778348646,0.00119368970687378,0.00674123722732964,0.000346164498632220,
                  -0.00120801661101190,-0.00868718303162277,0.000346164498632220,0.0130160227750794),4,4,byrow = TRUE)

## define prior
prior.num_params = 4;
prior.p1 = prior_mat[,1]
prior.p2 = prior_mat[,2]
prior.sampler<- function() apply(prior_mat,1,function(x) runif(1,x[1],x[2]))
prior.pdf <- function(theta) prod(exp(theta)/(1 + exp(theta))^2)
prior.trans_f<-function(theta) log((theta - prior.p1)/(prior.p2 - theta))
prior.trans_finv<-function(theta) (prior.p2*exp(theta) + prior.p1)/(1 + exp(theta))

theta = prior.trans_f(c(3,1,2,0.5))

if(exact) lval<-sum(dgk(y0, A=a, B=b, g=g, k=k,c=c,log=TRUE))

for(j in 1:5e4){
  
theta_1 = mvtnorm::rmvnorm(1,mean = theta,sigma = cov_rw)

prop = prior.trans_finv(theta_1)
a=prop[1]
b=prop[2]
g=prop[3]
k=prop[4]

if(exact)  {
  lval_1<-sum(dgk(y0, A=a, B=b, g=g, k=k,c=c,log=TRUE))
  logprior = log(prior.pdf(theta))
  logprior_1 = log(prior.pdf(theta_1))
  
  if(log(runif(1))<(lval_1-lval+logprior_1-logprior)){
    # if( abs(mean(xx)-0)<xx_q){
    lval=lval_1
    theta=theta_1
    cnt=cnt+1
  }
}else{
yy<-rgk(1e4, A=a, B=b, g=g, k=k,c=c)
# yy<-yeo.johnson(yy,lambda)
# yy<-(yy-muy)/sdy

clusterExport(cl,c(
  "yy"
),envir = environment())

f_est1_1<-unlist(parLapply(cl,y_ind,f_est))
if(length(f_est1_1)!=n) stop()
# f_est1_1<-unlist(parLapply(cl,1:length(y),f_est_orig))
(f_est1_0<-which(f_est1_1<0))
if(length(f_est1_0)<n/100){
  f_est1_1[f_est1_0]<-min(f_est1_1[f_est1_1>0])
  lval_1<-sum(log(f_est1_1))

# f_est1_1<-kde(x=yy,eval.points = y)$estimate
# (f_est1_0<-which(f_est1_1<0))
# if(length(f_est1_0)<n/100){
#   f_est1_1[f_est1_0]<-min(f_est1_1[f_est1_1>0])
#   lval_1<-sum(log(f_est1_1))

  
  if(!is.na(lval_1)){
    
    logprior = log(prior.pdf(theta))
    logprior_1 = log(prior.pdf(theta_1))
    
    if(log(runif(1))<(lval_1-lval+logprior_1-logprior)){
      # if( abs(mean(xx)-0)<xx_q){
      lval=lval_1
      theta=theta_1
      cnt=cnt+1
    }
  }
}
}

post_samp[j,]<-theta
if(j%%100==0) cat(paste0("iter=",j),"|",paste0("cnt=",cnt),"\n")
if(j%%10==0) print(round(c(theta),2))
}

if(!exact){
post_samp_f<-post_samp
for (j in 1:dim(post_samp_f)[1]) post_samp_f[j,] = prior.trans_finv(post_samp[j,])
}else{
post_samp_exact<-post_samp
for (j in 1:dim(post_samp_exact)[1]) post_samp_exact[j,] = prior.trans_finv(post_samp[j,])
}

stopCluster(cl)
}


### m/g/1



simulate_mg1<-function(T,theta){

  T<-T+1
  
w = rexp(T,theta[3])
s = runif(T,theta[1], theta[2])

x = numeric(T)
d = 0
v = 0;
for (i  in 1:T){
v = v + w[i]
d = d + s[i] + max(0,v-d)
x[i] = d
}

return(diff(x))
}

theta<-c(1,5,0.2)

n<-50
y0<-simulate_mg1(n,theta)
hist(y0)

muy<-mean(y0)
sdy<-sd(y0)
y<-(y0-muy)/sdy
y<-y0
hist(y)


q_y<-quantile(y,seq(0,1,by=0.1))
q_y<-quantile(y,c(0,0.1,0.25,0.75,0.9,1))
q_y<-quantile(y,c(0,1))
# q_y<-quantile(y,c(0,0.25,0.75,1))

yq_y<-NULL
for(i in 1:(length(q_y)-1)) yq_y[[i]]<-y[y>=q_y[i]&y<q_y[i+1]]
yq_y[[i]]<-c(yq_y[[i]],max(y))


lambda_y<-NULL
for(i in 1:length(yq_y)) lambda_y[[i]]<-yeojohnson(yq_y[[i]])$lambda

y_yj1_y<-
  y_yj1_mu<-
  y_yj1_sd<-NULL
for(i in 1:length(yq_y)){
  y_yj1_y[[i]]<-yeo.johnson(yq_y[[i]],lambda_y[[i]])
  y_yj1_mu[[i]]<-mean(y_yj1_y[[i]])
  y_yj1_sd[[i]]<-sd(y_yj1_y[[i]])
  y_yj1_y[[i]]<-(y_yj1_y[[i]]-y_yj1_mu[[i]])/y_yj1_sd[[i]]
} 

for(i in 1:length(y_yj1_y)) hist(y_yj1_y[[i]])

y_ind<-list()
y_len<-sapply(y_yj1_y,length)
cnt=1
for(i in 1:length(y_len)){  
  for(j in 1:(y_len[i])){
    y_ind[[cnt]]<-c(j,i)
    cnt=cnt+1
  }}

flambda<-function(y,lambda){
  if(y>=0) val<-abs((y+1)^(lambda-1))
  if(y<0) val<-abs((1-y)^(1-lambda))
  return(val)
}

R_y<-NULL
for(i in 1:length(yq_y)) R_y[[i]]<-floor(1+apply(as.matrix(yq_y[[i]]),1,
                                                 function(x) 10/
                                                   (y_yj1_sd[[i]]*flambda(x,lambda_y[[i]])) 
))



f_est<-function(arg1){
  
  i<-arg1[1]
  j<-arg1[2]
  yy_yj1_y<-yeo.johnson(yy,lambda_y[[j]])
  yy_yj1_y<-(yy_yj1_y-y_yj1_mu[[j]])/y_yj1_sd[[j]]
  
  f_est1<-fr_Rm_cpp(y_yj1_y[[j]][[i]],as.matrix(yy_yj1_y),
                    seq(5,15,length.out=10))
  
  return(median(f_est1))
  
}

f_est_orig<-function(arg1){
  
  i<-arg1[1]
  
  f_est1<-fr_Rm_cpp(y[i],as.matrix(yy),5:15)
  
  return(median(f_est1))
  
}


no.cores<-detectCores()-1
cl <- makeCluster(getOption("cl.cores", no.cores))

theta<-c(1,5,0.2)
yy<-simulate_mg1(1e4,theta)


clusterExport(cl,c("y","yq_y",
                   "y_yj1_y","yy",
                   "lambda_y","R_y",
                   "flambda","muy","sdy",
                   "y_yj1_mu","y_yj1_sd"
),
envir = environment()
)

# invisible(clusterEvalQ(cl, library(fourier)))
invisible(clusterEvalQ(cl, library(fourier)))
invisible(clusterEvalQ(cl, library(VGAM)))

f_est1_1<-unlist(parLapply(cl,y_ind,f_est))

(f_est1_0<-which(f_est1_1<0))
f_est1_1[f_est1_0]=min(f_est1_1[f_est1_1>0])
(lval<-sum(log(f_est1_1)))

# theta<-c(3,1,2,0.5)
# a=theta[1]
# b=theta[2]
# g=theta[3]
# k=theta[4]
# sum(dgk(y0, A=a, B=b, g=g, k=k,c=c,log = TRUE))
# lval<-sum(dgk(y0, A=a, B=b, g=g, k=k,c=c,log = TRUE))



post_samp<-matrix(0,1e4,3)
cnt=0

prior_mat<-matrix(0,3,2)
prior_mat[1,]<-c(0,min(y))
prior_mat[2,]<-c(0,10+min(y))
prior_mat[3,]<-c(0,0.5)

for(j in 1:1e4){
  
  theta_1<-rmvnorm(1,theta,c(.1^2,0.1^2,0.1^2)*diag(3))
  
  
  check_prior<-TRUE
  for(i in 1:3) {
    if((theta_1[i]<prior_mat[i,1]|theta_1[i]>prior_mat[i,2])) 
    {
      check_prior<-FALSE
      break
    }
  }
  
  
  if(check_prior){
    
    # for(i in 1:length(y_yj)) theta_list[[i]]<-c(y_yj[i],theta_1)
    # val1<-(unlist(parSapply(cl,theta_list,f_est)))
    # lval_1<-sum(log(val1))

    yy<-simulate_mg1(1e4,theta_1)
    # yy<-yeo.johnson(yy,lambda)
    # yy<-(yy-muy)/sdy

    clusterExport(cl,c(
      "yy"
    ),envir = environment())

    f_est1_1<-unlist(parLapply(cl,y_ind,f_est))
    if(length(f_est1_1)!=n) stop()
    # f_est1_1<-unlist(parLapply(cl,1:length(y),f_est_orig))
    (f_est1_0<-which(f_est1_1<0))
    if(length(f_est1_0)<n/100){
      f_est1_1[f_est1_0]<-min(f_est1_1[f_est1_1>0])
      lval_1<-sum(log(f_est1_1))
    
    
    if(!is.na(lval_1)){
      if(log(runif(1))<(lval_1-lval)){
        # if( abs(mean(xx)-0)<xx_q){
        lval=lval_1
        theta=theta_1
        cnt=cnt+1
      }
    }
  }
  }
  post_samp[j,]<-theta
  if(j%%100==0) cat(paste0("iter=",j),"|",paste0("cnt=",cnt),"\n")
  if(j%%10==0) print(round(c(theta),2))
}

stopCluster(cl)


### ricker 

## fourier

library(synlik)
library(matrixStats)
library(VGAM)
library(bestNormalize)
library(fourier)
library(parallel)
library(Rcpp)

set.seed(5)

n.t=100

# ricker_samp<-function(T,r,sigma,phi){
#   N<-numeric(T)
#   pop<-numeric(T)
#   N[1] = 7
#   for(t in 2:(T)){
#     e = rnorm(1,0,sigma)
#     if(t==1)   N[1] = 7
#     if(t>1) N[t] = r*N[t-1]*exp(-N[t-1]+e)
#     pop[t] = rpois(1,phi*N[t])
#     
#   }
#   return(list(pop=pop[1:(T)],N=N))
# }

cppFunction('List ricker_samp(int T, double r, double sigma, double phi) {

NumericVector N(T);
NumericVector pop(T);

N[0]=7.0;
pop[0]=R::rpois(phi*N[0]);
NumericVector e=Rcpp::rnorm(T,0,sigma);
for(int t = 1; t < T; ++t) {
N[t]=r*N[t-1]*exp(-N[t-1]+e[t]);
pop[t]=R::rpois(phi*N[t]);
}
    return List::create(_["pop"]=pop,_["N"]=N);
}')

#### Create synlik object
ricker_sl <- synlik(simulator = rickerSimul,
                    summaries = rickerStats,
                    param = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)),
                    extraArgs = list("nObs" = n.t, "nBurn" = 0, "initVal"=7.0,"randInit"=FALSE),
                    plotFun = function(input, ...){
                      plot(drop(input), type = 'l', ylab = "Pop", xlab = "Time", ...)
                    }
)

#### Simulate from the 
ricker_sl@data <- simulate(ricker_sl)
ricker_sl@extraArgs$obsData <- ricker_sl@data


y1<-c(ricker_sl@data)
# y1=ricker_samp(n.t,exp(3.8),0.3,10)$pop
y0<-y1+rnorm(length(y1),0,sd=.1)


muy<-mean(y0)
sdy<-sd(y0)
y<-(y0-muy)/sdy
y<-y0
hist(y)
plot.ts(y)
acf(y)
pacfy<-pacf(y)
cband<-qnorm((1 + 0.95)/2)/sqrt(length(y))
which.min((cband>abs(pacfy$acf))==0)-1



q_y<-quantile(y,seq(0,1,by=0.1))
q_y<-quantile(y,c(0,0.1,0.25,0.75,0.9,1))
q_y<-quantile(y,c(0,0.25,0.75,1))
q_y<-quantile(y,c(0,1))

yq_y<-NULL
for(i in 1:(length(q_y)-1)) yq_y[[i]]<-y[y>=q_y[i]&y<q_y[i+1]]
yq_y[[i]]<-c(yq_y[[i]],max(y))

yq_y_ind<-matrix(0,ncol = 2,n.t)
for(i in 1:n.t) {
  yq_y_ind_1<-which(sapply(yq_y,function(x) y[i] %in% x))
  yq_y_ind_2<-which(y[i]==yq_y[[yq_y_ind_1]])
  yq_y_ind[i,]<-c(yq_y_ind_1,yq_y_ind_2)
}


lambda_y<-NULL
for(i in 1:length(yq_y)) lambda_y[[i]]<-yeojohnson(yq_y[[i]])$lambda
# lambda_y[[1]]<-1

y_yj1_y<-
  y_yj1_mu<-
  y_yj1_sd<-NULL
for(i in 1:length(yq_y)){
  y_yj1_y[[i]]<-yeo.johnson(yq_y[[i]],lambda_y[[i]])
  y_yj1_mu[[i]]<-mean(y_yj1_y[[i]])
  y_yj1_sd[[i]]<-sd(y_yj1_y[[i]])
  # y_yj1_mu[[i]]<-0
  # y_yj1_sd[[i]]<-1
  y_yj1_y[[i]]<-(y_yj1_y[[i]]-y_yj1_mu[[i]])/y_yj1_sd[[i]]
}

for(i in 1:length(y_yj1_y)) hist(y_yj1_y[[i]])


y_ind_n<-y_ind_d<-list()
y_len<-sapply(y_yj1_y,length)

arx<-8

for(i in 1:(length(y)-arx)){  
  temp_ind<-NULL
  for(ii in i:(i+arx)){
    temp_ind<-rbind(temp_ind,c(yq_y_ind[ii,],ii))
  }
  y_ind_n[[i]]<-temp_ind
}


for(i in 1:(length(y)-(arx-1))){  
  temp_ind<-NULL
  for(ii in i:(i+(arx-1))){
    temp_ind<-rbind(temp_ind,c(yq_y_ind[ii,],ii))
  }
  y_ind_d[[i]]<-temp_ind
}


flambda<-function(y,lambda){
  if(y>=0) val<-abs((y+1)^(lambda-1))
  if(y<0) val<-abs((1-y)^(1-lambda))
  return(val)
}

R_y<-NULL
for(i in 1:length(yq_y)) R_y[[i]]<-floor(1+apply(as.matrix(yq_y[[i]]),1,
                                                 function(x) 10/
                                                   (y_yj1_sd[[i]]*flambda(x,lambda_y[[i]])) 
))



f_est<-function(arg1){
  
  n<-dim(arg1)[1]
  yy_yj1_y<-NULL
  y_yj1_y1<-NULL
  
  for(ii in 1:n){
    
    i<-arg1[ii,2]
    j<-arg1[ii,1]
    k<-arg1[ii,3]
    
    yy_yj1_y_temp<-yeo.johnson(yy[,k],lambda_y[[j]])
    yy_yj1_y_temp<-(yy_yj1_y_temp-y_yj1_mu[[j]])/y_yj1_sd[[j]]
    
    yy_yj1_y<-cbind(yy_yj1_y,yy_yj1_y_temp)
    y_yj1_y1<-c(y_yj1_y1,y_yj1_y[[j]][[i]])
    
  }
  
  
  if(n>8) f_est1<-fr_Rm_cpp(y_yj1_y1,as.matrix(yy_yj1_y),
                            seq(1,20,length.out=20))
  if(n<=8) f_est1<-fr_Rm_cpp(y_yj1_y1,as.matrix(yy_yj1_y),
                             seq(1,20,length.out=20))
  
  return((f_est1))
  
}

f_est_orig<-function(arg1){
  
  i<-arg1[1]
  
  f_est1<-fr_Rm_cpp(y[i],as.matrix(yy),5:15)
  
  return(median(f_est1))
  
}



no.cores<-detectCores()-1
# no.cores<-10
cl <- makeCluster(getOption("cl.cores", no.cores))


nsim<-1e4
theta <-  c(3.8,.3,10)

#### Create synlik object
ricker_sl1 <- synlik(simulator = rickerSimul,
                     summaries = rickerStats,
                     param = c(logR = 1, logSigma = 1, logPhi = 1),
                     extraArgs = list("nObs" = n.t, "nBurn" = 0, "initVal"=7.0,"randInit"=FALSE),
                     plotFun = function(input, ...){
                       plot(drop(input), type = 'l', ylab = "Pop", xlab = "Time", ...)
                     }
)

ricker_sl1@param<-c(logR = (th[1,1]), logSigma = th[2,1], logPhi = th[3,1])

yy=simulate(ricker_sl1, nsim)
# yy=t(replicate(ricker_samp(n.t,exp(th[1,1]),exp(th[2,1]),exp(th[3,1]))$pop,n=nsim))

yy=yy+rnorm(length(yy),0,sd = 0.1)

clusterExport(cl,c("y","yq_y",
                   "y_yj1_y","yy",
                   "lambda_y","R_y",
                   "flambda","muy","sdy",
                   "y_yj1_mu","y_yj1_sd"
),
envir = environment()
)

# invisible(clusterEvalQ(cl, library(fourier)))
invisible(clusterEvalQ(cl, library(fourier)))
invisible(clusterEvalQ(cl, library(VGAM)))

m<-10

temp_n_mat<-NULL
temp_n<-parLapply(cl,y_ind_n,f_est)
for(i in 1:length(y_ind_n)){
  temp_n_mat<-rbind(temp_n_mat,c(temp_n[[i]])/m)
  # plot.ts(temp_n[[i]])
}

temp_d_mat<-NULL
temp_d<-parLapply(cl,y_ind_d,f_est)
for(i in 1:length(y_ind_d)){
  temp_d_mat<-rbind(temp_d_mat,c(temp_d[[i]])/m)
  # plot.ts(temp_n[[i]])
}


for(jj in 1:m){
  temp_n<-parLapply(cl,y_ind_n,f_est)
  for(i in 1:length(y_ind_n)){
    temp_n_mat[i,]<-temp_n_mat[i,]+c(temp_n[[i]])/m
    # plot.ts(temp_n[[i]])
  }
  temp_d<-parLapply(cl,y_ind_d,f_est)
  for(i in 1:length(y_ind_d)){
    temp_d_mat[i,]<-temp_d_mat[i,]+c(temp_d[[i]])/m
    # plot.ts(temp_n[[i]])
  }
  
}


plot.ts(log(temp_n_mat[1,]),ylim=c(-20,3))
for(i in 1:length(y_ind_n)){
  lines(log(temp_n_mat[i,]))
}
which(is.nan(log(apply(temp_n_mat,1,function(x) median(x[3:6])))))
abline(h=log(apply(temp_n_mat,1,function(x) median(x[3:6]))),col=2)

plot.ts(log(temp_d_mat[1,]),ylim=c(-20,3),col=2)
for(i in 1:length(y_ind_n)){
  lines(log(temp_d_mat[i,]),col=2)
}
which(is.nan(log(apply(temp_d_mat,1,function(x) median(x[3:7])))))
abline(h=log(apply(temp_d_mat,1,function(x) median(x[3:7]))),col=3)


y_ind_n_r<-apply(temp_n_mat,1,function(x) max(3,min(5,which(x[1:5]<0))))
y_ind_d_r<-apply(temp_d_mat,1,function(x) max(3,min(5,which(x[1:5]<0))))

y_ind_n<-y_ind_d<-list()
for(i in 1:(length(y)-arx)){  
  temp_ind<-NULL
  for(ii in i:(i+arx)){
    temp_ind<-rbind(temp_ind,c(yq_y_ind[ii,],ii,y_ind_n_r[i]))
  }
  y_ind_n[[i]]<-temp_ind
}


for(i in 1:(length(y)-(arx-1))){  
  temp_ind<-NULL
  for(ii in i:(i+(arx-1))){
    temp_ind<-rbind(temp_ind,c(yq_y_ind[ii,],ii,y_ind_d_r[i]))
  }
  y_ind_d[[i]]<-temp_ind
}

stopCluster(cl)


f_est<-function(arg1){
  
  n<-dim(arg1)[1]
  yy_yj1_y<-NULL
  y_yj1_y1<-NULL
  
  for(ii in 1:n){
    
    i<-arg1[ii,2]
    j<-arg1[ii,1]
    k<-arg1[ii,3]
    
    yy_yj1_y_temp<-yeo.johnson(yy[,k],lambda_y[[j]])
    yy_yj1_y_temp<-(yy_yj1_y_temp-y_yj1_mu[[j]])/y_yj1_sd[[j]]
    
    yy_yj1_y<-cbind(yy_yj1_y,yy_yj1_y_temp)
    y_yj1_y1<-c(y_yj1_y1,y_yj1_y[[j]][[i]])
    
  }
  
  
  # r1=min(3,mean(arg1[,4]))
  r2=mean(arg1[,4])
  
  if(n>8) f_est1<-fr_Rm_cpp(y_yj1_y1,as.matrix(yy_yj1_y),
                            seq(3,6,length.out=5))
  if(n<=8) f_est1<-fr_Rm_cpp(y_yj1_y1,as.matrix(yy_yj1_y),
                             seq(3,6,length.out=5))
  
  return(median(f_est1))
  
}


no.cores<-detectCores()-1
# no.cores<-10
cl <- makeCluster(getOption("cl.cores", no.cores))


m<-10
NN<-5e3

nsim<-1e4
theta <-  c(3.8,.3,10)
th <- matrix(0,3,NN) ## storage for chain results
th[,1] <- c(theta[1],log(theta[2]),log(theta[3]))        ## initial state
lv<-numeric(NN)

propCov = diag(c(0.1, 0.1, 0.1))^2
cholFact <- t(chol(unname(propCov)))
nPar <- 3

#### Create synlik object
ricker_sl1 <- synlik(simulator = rickerSimul,
                     summaries = rickerStats,
                     param = c(logR = 1, logSigma = 1, logPhi = 1),
                     extraArgs = list("nObs" = n.t, "nBurn" = 0, "initVal"=7.0,"randInit"=FALSE),
                     plotFun = function(input, ...){
                       plot(drop(input), type = 'l', ylab = "Pop", xlab = "Time", ...)
                     }
)

ricker_sl1@param<-c(logR = (th[1,1]), logSigma = th[2,1], logPhi = th[3,1])

f_est1_n1<-f_est1_d1<-numeric(n.t-(arx-1))
for(jj in 1:m){
  
  yy=simulate(ricker_sl1, nsim)
  # yy=t(replicate(ricker_samp(n.t,exp(th[1,1]),exp(th[2,1]),exp(th[3,1]))$pop,n=nsim))
  
  yy=yy+rnorm(length(yy),0,sd = 0.1)
  
  clusterExport(cl,c("y","yq_y",
                     "y_yj1_y","yy",
                     "lambda_y","R_y",
                     "flambda","muy","sdy",
                     "y_yj1_mu","y_yj1_sd"
  ),
  envir = environment()
  )
  
  invisible(clusterEvalQ(cl, library(fourier)))
  invisible(clusterEvalQ(cl, library(VGAM)))
  
  
  f_est1_n1_temp<-unlist(parLapply(cl,y_ind_n,f_est))
  f_est1_d1_temp<-unlist(parLapply(cl,y_ind_d,f_est))
  if(length(f_est1_n1_temp)!=(n.t-arx)) stop()
  if(length(f_est1_d1_temp)!=(n.t-(arx-1))) stop()
  
  f_est1_n1<-f_est1_n1+c(f_est(y_ind_d[[1]]),f_est1_n1_temp)/m
  f_est1_d1<-f_est1_d1+f_est1_d1_temp/m
}

(f_est1_n0<-which(f_est1_n1<0))
(f_est1_d0<-which(f_est1_d1<0))

f_est1_n1[f_est1_n0]=min(f_est1_n1[f_est1_n1>0])
f_est1_d1[f_est1_d0]=min(f_est1_d1[f_est1_d1>0])

(lval<-sum(log(f_est1_n1))-sum(log(f_est1_d1)))
lv[1]<-lval


cnt=0
cntx<-0
cntxx<-0
check_cntxx<-0

for(j in 2:NN){
  
  if(cntx>5){
    ricker_sl1@param<-c(logR = (th[1,j-1]), logSigma = th[2,j-1], logPhi = th[3,j-1])
    cntxx<-1
    f_est1_n1_total<-f_est1_d1_total<-numeric(n.t-(arx-1))
    while(cntxx>0){
      if(cntxx<=3){
        f_est1_n1<-f_est1_d1<-numeric(n.t-(arx-1))
        for(jj in 1:m){
          
          yy=simulate(ricker_sl1, nsim)
          # yy=t(replicate(ricker_samp(n.t,exp(th[1,j-1]),exp(th[2,j-1]),exp(th[3,j-1]))$pop,n=nsim))
          
          yy=yy+rnorm(length(yy),0,sd = 0.1)
          
          clusterExport(cl,c("yy"
          ),
          envir = environment()
          )
          
          invisible(clusterEvalQ(cl, library(fourier)))
          invisible(clusterEvalQ(cl, library(VGAM)))
          
          
          f_est1_n1_temp<-unlist(parLapply(cl,y_ind_n,f_est))
          f_est1_d1_temp<-unlist(parLapply(cl,y_ind_d,f_est))
          
          if(length(f_est1_n1_temp)!=(n.t-arx)) stop()
          if(length(f_est1_d1_temp)!=(n.t-(arx-1))) stop()
          
          f_est1_n1<-f_est1_n1+c(f_est(y_ind_d[[1]]),f_est1_n1_temp)
          f_est1_d1<-f_est1_d1+f_est1_d1_temp
        }
        
        f_est1_n1_total=(f_est1_n1_total*(m*(cntxx-1))+f_est1_n1)/(m*cntxx)
        f_est1_d1_total=(f_est1_d1_total*(m*(cntxx-1))+f_est1_d1)/(m*cntxx)
        # f_est1_n1_total=f_est1_n1/m
        # f_est1_d1_total=f_est1_d1/m
        
        # f_est1_1<-unlist(parLapply(cl,1:length(y),f_est_orig))
        (f_est1_n0<-which(f_est1_n1_total<=0))
        (f_est1_d0<-which(f_est1_d1_total<=0))
        
        if(length(f_est1_n0)<2&length(f_est1_d0)<2){
          f_est1_n1_total[f_est1_n0]=min(f_est1_n1_total[f_est1_n1_total>0])
          f_est1_d1_total[f_est1_d0]=min(f_est1_d1_total[f_est1_d1_total>0])
          
          lval<-sum(log(f_est1_n1_total))-sum(log(f_est1_d1_total))
          cntx=0
          cntxx=0
        }else{
          cntxx<-cntxx+1
        }
      }else{
        lval=min(lv[1:(j-1)])
        cntx=0
        cntxx=0
        check_cntxx=check_cntxx+1
        print(check_cntxx)
      }
    }
  }
  
  
  th[,j] <- th[,j-1]
  
  # th_ind<-sample(1:3,1)
  # # th_ind<-2
  # pert <- rnorm(1)
  # th[th_ind,j] <- th[th_ind,j] + as.vector(cholFact[th_ind,th_ind] %*% pert)
  
  pert <- rnorm(nPar)
  th[,j] <- th[,j] + as.vector(cholFact %*% pert)
  
  # if(theta_1[1]>0){
  
  ricker_sl1@param<-c(logR = (th[1,j]), logSigma = th[2,j], logPhi = th[3,j])
  
  f_est1_n1<-f_est1_d1<-numeric(n.t-(arx-1))
  for(jj in 1:m){
    
    yy=simulate(ricker_sl1, nsim)
    # yy=t(replicate(ricker_samp(n.t,exp(th[1,j]),exp(th[2,j]),exp(th[3,j]))$pop,n=nsim))
    
    yy=yy+rnorm(length(yy),0,sd = 0.1)
    
    clusterExport(cl,c("yy"
    ),
    envir = environment()
    )
    
    invisible(clusterEvalQ(cl, library(fourier)))
    invisible(clusterEvalQ(cl, library(VGAM)))
    
    
    f_est1_n1_temp<-unlist(parLapply(cl,y_ind_n,f_est))
    f_est1_d1_temp<-unlist(parLapply(cl,y_ind_d,f_est))
    if(length(f_est1_n1_temp)!=(n.t-arx)) stop()
    if(length(f_est1_d1_temp)!=(n.t-(arx-1))) stop()
    
    f_est1_n1<-f_est1_n1+c(f_est(y_ind_d[[1]]),f_est1_n1_temp)/m
    f_est1_d1<-f_est1_d1+f_est1_d1_temp/m
  }
  
  (f_est1_n0<-which(f_est1_n1<=0))
  (f_est1_d0<-which(f_est1_d1<=0))
  
  if(length(f_est1_n0)<2&&length(f_est1_d0)<2){
    f_est1_n1[f_est1_n0]=min(f_est1_n1[f_est1_n1>0])
    f_est1_d1[f_est1_d0]=min(f_est1_d1[f_est1_d1>0])
    
    lval_1<-sum(log(f_est1_n1))-sum(log(f_est1_d1))
    
    
    
    if(log(runif(1))<(lval_1-lval)){
      # if( abs(mean(xx)-0)<xx_q){
      lval=lval_1
      # theta=theta_1
      cnt=cnt+1
      cntx=0
      cntxx=0
    }else{
      th[,j] <- th[,j-1]
      cntx=cntx+1 
    }
  }else{
    th[,j] <- th[,j-1]
  }
  # }
  
  lv[j]<-lval
  
  # if(j%%10==0){
  cat(paste0("iter=",j),"|",paste0("theta=",c(round(th[1,j],3),round(exp(th[2:3,j]),3))), "|",   
      paste0("accept=",round(cnt/j,3)),"\n")
  # }
}

stopCluster(cl)



## sl

n.t<-100
#### Create synlik object
ricker_sl <- synlik(simulator = rickerSimul,
                    summaries = rickerStats,
                    param = c(logR = 3.8, logSigma = log(0.3), logPhi = log(10)),
                    extraArgs = list("nObs" = n.t, "nBurn" = 0,"initVal"=7.0,"randInit" = FALSE),
                    plotFun = function(input, ...){
                      plot(drop(input), type = 'l', ylab = "Pop", xlab = "Time", ...)
                    }
)

#### Simulate from the object
ricker_sl@data <- simulate(ricker_sl)
ricker_sl@extraArgs$obsData <- ricker_sl@data

NN<-5e3
# samp_sl<-function(NN){
nsim=1e4
n.t=n.t
nn<-10
burn=0
theta <-  c(3.8,.3,10)
# theta[1]<-runif(1,3,4)
# theta[2]<-runif(1,.1,.4)
# theta[3]<-runif(1,7,12)

propCov = diag(c(0.1, 0.2, 0.1))^2
cholFact <- t(chol(unname(propCov)))
nPar <- 3



llr <- rep(NA,NN)         ## storage for log synthetic likelihood
th <- matrix(0,3,NN) ## storage for chain results
se <- rep(NA,NN)  
th[,1] <- c(theta[1],log(theta[2]),log(theta[3]))        ## initial state


library(parallel)


cl <- makeCluster(getOption("cl.cores", detectCores()-1))

clusterExport(cl,c("ricker_sl","th","nsim"))
invisible(clusterEvalQ(cl, library(fourier)))
invisible(clusterEvalQ(cl, library(synlik)))


# sl=0
#   for(ii in 1:nn){
#         sl=sl+slik(ricker_sl, 
#          param = c(logR = (th[1,1]), logSigma = th[2,1], logPhi = th[3,1]),
#          nsim   = nsim)/nn
#     
#   }
# val0<-sl

sl_ll<-function(th){
  slik(ricker_sl, 
       param = c(logR = (th[1]), logSigma = th[2], logPhi = th[3]),
       nsim   = nsim)
}


val0<-mean(parSapply(cl, replicate(nn,list((th[,1]))),sl_ll))

(llr[1] <-(val0))

prop.sd <- c(.02,.1,.05)  #on log scale
reject <- 0
uni <- runif(NN)

cntx<-0

for (i in 2:(NN+burn)) { ## MCMC loop
  
  
  ##prev
  if(cntx>5){
    # sl=0
    #   for(ii in 1:nn){
    #     sl=sl+slik(ricker_sl,
    #                param = c(logR = (th[1,i-1]), logSigma = th[2,i-1], logPhi = th[3,i-1]),
    #                nsim   = nsim)/nn
    # 
    #   }
    # 
    #   val0<-sl
    
    val0<-mean(parSapply(cl, replicate(nn,list((th[,i-1]))),sl_ll))
    
    cntx<-0
  }
  
  
  th[,i] <- th[,i-1]
  
  # th_ind<-sample(1:3,1)
  # # th_ind<-2
  # pert <- rnorm(1)
  # th[th_ind,i] <- th[th_ind,i] + as.vector(cholFact[th_ind,th_ind] %*% pert)
  
  pert <- rnorm(nPar)
  th[,i] <- th[,i] + as.vector(cholFact %*% pert)
  
  
  # sl=0
  # for(ii in 1:nn){
  #   sl=sl+slik(ricker_sl,
  #              param = c(logR = (th[1,i]), logSigma = th[2,i], logPhi = th[3,i]),
  #              nsim   = nsim)/nn
  # 
  # }
  # 
  # val1<-sl
  
  
  val1<-mean(parSapply(cl, replicate(nn,list((th[,i]))),sl_ll))
  
  
  alpha <- min(1,exp(val1-val0))
  
  
  if (uni[i]>alpha) { ## reject
    th[,i] <- th[,i-1]
    llr[i] <- val0
    reject <- reject + 1
    cntx=cntx+1
  }else{
    cntx=0
    llr[i] <- val1
  }
  
  
  
  
  cat(paste0("iter=",i),"|",paste0("theta=",c(round(th[1,i],3),round(exp(th[2:3,i]),3))), "|",   
      paste0("accept=",round(1-reject/i,3)),"\n")
  
  
}

#   if(j%%1e2==0) cat(paste0("iter=",j),"|",paste0("cnt=",cnt),"\n")
etime<-Sys.time() 
#   return(list(theta=cbind(th[1,],exp(th[2,]),exp(th[3,])),llr=llr))
# }

stopCluster(cl)

## particle

library(pomp)
library(dplyr)

dat<-data.frame(cbind(times=1:n.t,pop=1:n.t))


ricker_pf <- pomp(dat,times="times",t0=1)
plot(ricker_pf)

stochStep <- Csnippet("
  e = rnorm(0,sigma);
  N = r*N*exp(-c*N+e);
")
pomp(ricker_pf,rprocess=discrete_time(step.fun=stochStep,delta.t=1),
     paramnames=c("r","c","sigma"),statenames=c("N","e")) -> ricker_pf

rmeas <- Csnippet("pop = rpois(phi*N);")
dmeas <- Csnippet("lik = dpois(pop,phi*N,give_log);")
pomp(ricker_pf,rmeasure=rmeas,dmeasure=dmeas,statenames=c("N"),paramnames=c("phi")) -> ricker_pf

skel <- Csnippet("DN = r*N*exp(-c*N);")
ricker_pf <- pomp(ricker_pf,skeleton=map(skel),paramnames=c("r","c"),statenames=c("N"))

coef(ricker_pf) <- c(N_0=1,e_0=0,r=exp(3.8),c=1,sigma=0.3,phi=10)

{
  ricker_pf %>%
    simulate(
      rinit=Csnippet("
        e = rnorm(0,sigma);
            // N = r*N_0*exp(-c*N_0+e);
             N = 7.0;
      for (int i=0; i < 100; i++) {
        e = rnorm(0,sigma);
  N = r*N*exp(-c*N+e);
      }
    "),
      paramnames=c("N_0","e_0","sigma","phi","r","c"),
      statenames=c("N","e")
    ) -> ricker_pf
  
}

x <- simulate(ricker_pf)
plot(x)
# x@data
# ricker_pf@data<-x@data

rname<-rownames(ricker_pf@data)
ricker_pf@data<-(y)
rownames(ricker_pf@data)<-rname



# yy <- simulate(ricker_pf,nsim = 1e3)
# dim(sapply(yy, function(x) x@data))


# prior_mat<-matrix(0,3,2)
# prior_mat[1,]<-c(2,5)
# prior_mat[2,]<-c(-3,-0.22)
# prior_mat[3,]<-c(1.61,3)

prior_mat<-matrix(0,3,2)
prior_mat[1,]<-c(1,10)
prior_mat[2,]<-log(c(0.0001,2))
prior_mat[3,]<-log(c(1,30))

NN<-5000

  uni<-runif(NN)
  nsim=1e4
  n.t=n.t
  nn<-10
  burn=0
  theta <-  c(apply(prior_mat,1,mean)[1],c(exp(apply(prior_mat,1,mean)[2:3])))
  # theta[1]<-runif(1,3,4)
  # theta[2]<-runif(1,.1,.4)
  # theta[3]<-runif(1,7,9)
  theta<-c(3.8,0.3,10.0)
  
  propCov = diag(c(0.1, 0.1, 0.1))^2
  cholFact <- t(chol(unname(propCov)))
  nPar <- 3
  
  
  
  llr <- rep(NA,NN)         ## storage for log synthetic likelihood
  th <- matrix(0,3,NN) ## storage for chain results
  se <- rep(NA,NN)  
  th[,1] <- c(theta[1],log(theta[2]),log(theta[3]))        ## initial state
  
  
  coef(ricker_pf) <- c(N_0=1,e_0=0,r=exp(th[1,1]),c=1,sigma=exp(th[2,1]),phi=exp(th[3,1]))
  plink<-0
  for(ii in 1:nn){
    pf <- pfilter(ricker_pf,Np=nsim)	## use 1000 particles
    plink=plink+logLik(pf)/nn
  }
  
  pf <- pfilter(ricker(r=exp(th[1,1]),sigma = exp(th[2,1]),phi = exp(th[3,1])),Np=nsim)
  

  llr[1]<-plink
  
  reject=0
  
  for (i in 2:(NN+burn)) { ## MCMC loop
    
    ## prev
    
    coef(ricker_pf) <- c(N_0=1,e_0=0,r=exp(th[1,i-1]),c=1,sigma=exp(th[2,i-1]),phi=exp(th[3,i-1]))
    plink<-0
    for(ii in 1:nn){
      pf <- pfilter(ricker_pf,Np=nsim)	
      plink=plink+logLik(pf)/nn
    }
    
    llr[i-1]<-plink
    
    ## new
    
    th[,i] <- th[,i-1]
    
    # th_ind<-sample(1:3,1)
    # # th_ind<-2
    # pert <- rnorm(1)
    # th_temp<-th[th_ind,i-1] + as.vector(cholFact[th_ind,th_ind] %*% pert)
    # if(th_temp>prior_mat[th_ind,1]&&th_temp<prior_mat[th_ind,2]) {
    #
    #   th[th_ind,i] <- th_temp
    
    pert <- rnorm(nPar)
    th_temp <- th[,i] + as.vector(cholFact %*% pert)
    
    if(
      (th_temp[1]>prior_mat[1,1]&&th_temp[1]<prior_mat[1,2])&
      (th_temp[2]>prior_mat[2,1]&&th_temp[2]<prior_mat[2,2])&
      (th_temp[3]>prior_mat[3,1]&&th_temp[3]<prior_mat[3,2])
    ) {
      
      th[,i] <- th_temp
      
      
      # pert <- rnorm(nPar)
      # th[,i] <- th[,i] + as.vector(cholFact %*% pert)
      
      coef(ricker_pf) <- c(N_0=1,e_0=0,r=exp(th[1,i]),c=1,sigma=exp(th[2,i]),phi=exp(th[3,i]))
      plink<-0
      for(ii in 1:nn){
        pf <- pfilter(ricker_pf,Np=nsim)	
        plink=plink+logLik(pf)/nn
      }
      
      
      llr[i]<-plink
      
      
      
      RR<-(llr[i]-llr[i-1])
      
      alpha <- min(1,exp(RR))
      
      if (uni[i]>alpha) { ## reject
        th[,i] <- th[,i-1]
        llr[i] <- llr[i-1]
        reject <- reject + 1
      }
    }else{ ## reject
      th[,i] <- th[,i-1]
      llr[i] <- llr[i-1]
      reject <- reject + 1
      
    }
    
    
    cat(paste0("iter=",i),"|",paste0("theta=",c(round(th[1,i],3),round(exp(th[2:3,i]),3))), "|",   
        paste0("accept=",round(1-reject/i,3)),"\n")
    
    
  }








