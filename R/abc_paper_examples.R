library(fourier)
library(ks)
library(parallel)

############################################
############################################
### univariate mixture
############################################
############################################

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


############################################
############################################
### multivariate mixture
############################################
############################################

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


############################################
############################################
### g-and-k
############################################
############################################

library(gk)
library(mvtnorm)
library(bestNormalize)
library(VGAM)
library(fourier)
library(parallel)

### data needs to be imported from https://github.com/cdrovandi/ABC-dist-compare
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
  cl <- makeCluster(getOption("cl.cores", no.cores))
  
  theta<-c(3,1,2,0.5)
  a=theta[1]
  b=theta[2]
  g=theta[3]
  k=theta[4]
  c=0.8
  yy<-rgk(1e4, A=a, B=b, g=g, k=k,c=c)
  # yy<-yeo.johnson(yy,lambda)

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


############################################
### multivariate g-and-k
############################################


library(gk)
library(mvtnorm)

set.seed(1)
hyper <- runif(1, -10, 10)
As <- rnorm(50, hyper, 1)
B <- runif(1, 0, 1)
g <- runif(1, 0, 1)
k <- runif(1, 0, 1)

X <- matrix(NA, ncol = 50, nrow = 20)
for (i in 1:50) {
  X[, i] <- rgk(20, As[i], B, g, k)
  hist(X[,i])
}


###fourier

library(fourier)
library(bestNormalize)
library(VGAM)
library(parallel)


y_yj1_y_list<-
  y_yj1_mu_list<-
  y_yj1_sd_list<-
  y_ind_list<-
  lambda_y_list<-
  list()

for(kk in 1:50){
  
  y<-X[,kk]
  q_y<-quantile(y,c(0,0.25,0.75,1))
  
  yq_y<-NULL
  for(i in 1:(length(q_y)-1)) yq_y[[i]]<-y[y>=q_y[i]&y<q_y[i+1]]
  yq_y[[i]]<-c(yq_y[[i]],max(y))
  
  lambda_y<-NULL
  for(i in 1:length(yq_y)) lambda_y[[i]]<-yeojohnson(yq_y[[i]])$lambda
  lambda_y_list[[kk]]<-lambda_y
  
  y_yj1_y<-
    y_yj1_mu<-
    y_yj1_sd<-NULL
  for(i in 1:length(yq_y)){
    y_yj1_y[[i]]<-yeo.johnson(yq_y[[i]],lambda_y[[i]])
    y_yj1_mu[[i]]<-mean(y_yj1_y[[i]])
    y_yj1_sd[[i]]<-sd(y_yj1_y[[i]])
    y_yj1_y[[i]]<-(y_yj1_y[[i]]-y_yj1_mu[[i]])/y_yj1_sd[[i]]
  } 
  
  y_yj1_y_list[[kk]]<-y_yj1_y
  y_yj1_mu_list[[kk]]<-y_yj1_mu
  y_yj1_sd_list[[kk]]<-y_yj1_sd
  
  
  y_ind<-list()
  y_len<-sapply(y_yj1_y,length)
  cnt=1
  for(i in 1:length(y_len)){  
    for(j in 1:(y_len[i])){
      y_ind[[cnt]]<-c(j,i)
      cnt=cnt+1
    }}
  
  y_ind_list[[kk]]<-y_ind
  
}




f_est<-function(arg1){
  
  i<-arg1[1]
  j<-arg1[2]
  yy_yj1_y<-yeo.johnson(yy,lambda_y[[j]])
  yy_yj1_y<-(yy_yj1_y-y_yj1_mu[[j]])/y_yj1_sd[[j]]
  
  f_est1<-fr_Rm_cpp(y_yj1_y[[j]][[i]],as.matrix(yy_yj1_y),
                    seq(5,15,length.out=10))
  
  return(median(f_est1))
  
}


## Use option cl.cores to choose an appropriate cluster size.
no.cores<-detectCores()-1
cl <- makeCluster(getOption("cl.cores", no.cores))


set.seed(1)
m<-1e4
mt<-1e3
mm<-10
theta_post<-matrix(0,nrow = m,50)
alpha_post<-numeric(m)
fval_post<-numeric(m)

alpha<-mean(X)
theta<-apply(X,2,mean)

f_est1_1<-NULL

for(kk in 1:50){
  y_yj1_y<-y_yj1_y_list[[kk]]
  lambda_y<-lambda_y_list[[kk]]
  y_yj1_mu<-y_yj1_mu_list[[kk]]
  y_yj1_sd<-y_yj1_sd_list[[kk]]
  yy<-rgk(mt, theta[kk], B, g, k)
  
  clusterExport(cl,c("y_yj1_y","yy",
                     "lambda_y",
                     "y_yj1_mu","y_yj1_sd"
  ),
  envir = environment()
  )
  
  invisible(clusterEvalQ(cl, library(fourier)))
  invisible(clusterEvalQ(cl, library(VGAM)))
  invisible(clusterEvalQ(cl, library(gk)))
  
  
  f_est1_1<-c(f_est1_1,unlist(parLapply(cl,y_ind_list[[kk]],f_est)))
  
}

f_est1_1_0<-which(f_est1_1<0)
length(f_est1_1_0)

if(length(f_est1_1_0)<10){
  f_est1_1[f_est1_1_0]<-min(f_est1_1[f_est1_1>0])
}

(lval<-sum(log(f_est1_1)))
fval<-lval+sum(dnorm(theta, alpha,log=TRUE))

cnt=0
cntx=0
cntxx=2

for(i in 1:m){
  
  ### old
  
  if(cntx>5){
    
    alpha_1<-alpha
    theta_1<-theta
    
    f_est1_1<-NULL
    
    for(kk in 1:50){
      y_yj1_y<-y_yj1_y_list[[kk]]
      lambda_y<-lambda_y_list[[kk]]
      y_yj1_mu<-y_yj1_mu_list[[kk]]
      y_yj1_sd<-y_yj1_sd_list[[kk]]
      yy<-rgk(mt, theta_1[kk], B, g, k)
      
      clusterExport(cl,c("y_yj1_y","yy",
                         "lambda_y",
                         "y_yj1_mu","y_yj1_sd"
      ),
      envir = environment()
      )
      
      invisible(clusterEvalQ(cl, library(fourier)))
      invisible(clusterEvalQ(cl, library(VGAM)))
      invisible(clusterEvalQ(cl, library(gk)))
      
      
      f_est1_1<-c(f_est1_1,unlist(parLapply(cl,y_ind_list[[kk]],f_est)))
    }
    
    f_est1_1_0<-which(f_est1_1<0)
    if(length(f_est1_1_0)<10){
      f_est1_1[f_est1_1_0]<-min(f_est1_1[f_est1_1>0])
      (lval<-sum(log(f_est1_1)))
      fval<-lval+sum(dnorm(theta_1, alpha_1,log=TRUE))
      cntxx=2
    }else{
      fval<-fval_post[i-1]-log(cntxx)
      cntxx=cntxx+1
    }
    cntx=0
  }
  
  ### new  
  
  alpha_1<-rnorm(1,alpha,0.05)
  theta_1<-rnorm(50,theta,0.05)
  
  f_est1_1<-NULL
  
  for(kk in 1:50){
    y_yj1_y<-y_yj1_y_list[[kk]]
    lambda_y<-lambda_y_list[[kk]]
    y_yj1_mu<-y_yj1_mu_list[[kk]]
    y_yj1_sd<-y_yj1_sd_list[[kk]]
    yy<-rgk(mt, theta_1[kk], B, g, k)
    
    clusterExport(cl,c("y_yj1_y","yy",
                       "lambda_y",
                       "y_yj1_mu","y_yj1_sd"
    ),
    envir = environment()
    )
    
    # invisible(clusterEvalQ(cl, library(mypackage)))
    invisible(clusterEvalQ(cl, library(fourier)))
    invisible(clusterEvalQ(cl, library(VGAM)))
    invisible(clusterEvalQ(cl, library(gk)))
    
    
    f_est1_1<-c(f_est1_1,unlist(parLapply(cl,y_ind_list[[kk]],f_est)))
    
  }
  
  f_est1_1_0<-which(f_est1_1<0)
  if(length(f_est1_1_0)<10){
    f_est1_1[f_est1_1_0]<-min(f_est1_1[f_est1_1>0])
    (lval<-sum(log(f_est1_1)))
    fval1<-lval+sum(dnorm(theta_1, alpha_1,log=TRUE))
    
    if(log(runif(1))<(fval1-fval)){
      alpha=alpha_1
      theta=theta_1
      fval<-fval1
      cnt=cnt+1
      cntx=0
      cntxx=2
    }else{
      cntx=cntx+1
    }
  }
  
  alpha_post[i]<-alpha
  theta_post[i,]<-theta
  fval_post[i]<-fval
  
  print(paste0("iter=",i,"; alpha=",round(alpha,3),"; accept=",round(cnt/i,2)))
  
}

stopCluster(cl)


### exact

set.seed(5)
m<-5e4

theta_post_true<-matrix(0,nrow = m,50)
alpha_post_true<-numeric(m)
fval_post_true<-numeric(m)

alpha<-mean(X)
theta<-apply(X,2,mean)

lval=0
for(ii in 1:50){
  lval=lval+sum(dgk(X[,ii],theta[ii],B,g,k,log = TRUE))
}

fval<-lval+sum(dnorm(theta, alpha,log=TRUE))


cntt=0

for(i in 1:m){
  
  
  alpha_1<-rnorm(1,alpha,0.02)
  theta_1<-rnorm(50,theta,0.02)
  
  lval=0
  for(ii in 1:50){
    lval=lval+sum(dgk(X[,ii],theta_1[ii],B,g,k,log = TRUE))
  }
  
  fval1<-lval+sum(dnorm(theta_1, alpha_1,log=TRUE))
  
  if(log(runif(1))<(fval1-fval)){
    alpha=alpha_1
    theta=theta_1
    fval<-fval1
    cntt=cntt+1
  }
  
  
  alpha_post_true[i]<-alpha
  theta_post_true[i,]<-theta
  fval_post_true[i]<-fval
  
  print(paste0("iter=",i,"; alpha=",round(alpha,3),"; accept=",round(cntt/i,2)))
  
}


############################################
############################################
### m/g/1
############################################
############################################



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

    
    yy<-simulate_mg1(1e4,theta_1)

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



############################################
############################################
### ricker 
############################################
############################################

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

#### Simulate from the model
ricker_sl@data <- simulate(ricker_sl)
ricker_sl@extraArgs$obsData <- ricker_sl@data


y1<-c(ricker_sl@data)
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
                            seq(3,6,length.out=10))
  if(n<=8) f_est1<-fr_Rm_cpp(y_yj1_y1,as.matrix(yy_yj1_y),
                             seq(3,6,length.out=10))
  
  return(median(f_est1))
  
}


no.cores<-detectCores()-1
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
cntxx=2

for(j in 2:NN){
  
  if(cntx>5){
    ricker_sl1@param<-c(logR = (th[1,j-1]), logSigma = th[2,j-1], logPhi = th[3,j-1])
    
    f_est1_n1<-f_est1_d1<-numeric(n.t-(arx-1))
    for(jj in 1:m){
      
      yy=simulate(ricker_sl1, nsim)
      
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
      
      lval<-sum(log(f_est1_n1))-sum(log(f_est1_d1))
      cntxx=2
    }else{
      lval=lv[j-1]-log(cntxx)
      cntxx=cntxx+1
    }
    
  }
  
  
  th[,j] <- th[,j-1]
  
  # th_ind<-sample(1:3,1)
  # # th_ind<-2
  # pert <- rnorm(1)
  # th[th_ind,j] <- th[th_ind,j] + as.vector(cholFact[th_ind,th_ind] %*% pert)
  
  pert <- rnorm(nPar)
  th[,j] <- th[,j] + as.vector(cholFact %*% pert)
  
  
  ricker_sl1@param<-c(logR = (th[1,j]), logSigma = th[2,j], logPhi = th[3,j])
  
  f_est1_n1<-f_est1_d1<-numeric(n.t-(arx-1))
  for(jj in 1:m){
    
    yy=simulate(ricker_sl1, nsim)

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
      lval=lval_1
      cnt=cnt+1
      cntx=0
      cntxx=2
    }else{
      th[,j] <- th[,j-1]
      cntx=cntx+1 
    }
  }else{
    th[,j] <- th[,j-1]
  }
  # }
  
  lv[j]<-lval
  
  cat(paste0("iter=",j),"|",paste0("theta=",c(round(th[1,j],3),round(exp(th[2:3,j]),3))), "|",   
      paste0("accept=",round(cnt/j,3)),"\n")
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
nsim=1e4
nn<-10
burn=0
theta <-  c(3.8,.3,10)


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

stopCluster(cl)
