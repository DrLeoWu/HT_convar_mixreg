##################################################################################################################
##################################################################################################################
####################   Simulation investigating the hypothesis testing on coefficients of mixture of regressions 
##################################################################################################################
##################################################################################################################
library(mixtools)
library(sn)
library(matrixcalc)
library(Matrix)
library(ks)
library(pracma)
library(interp)
library(parallel)

emopt_pi <- function(initl,X,w){ # initl=initial value for beta, X is n by p, w is n by 1
  
  emobj_pi <- function(beta_pi){
    nsamp <- nrow(X)
    return(-c(t(w)%*%X%*%beta_pi)/nsamp + sum(log(1+exp(X%*%beta_pi)))/nsamp)
  }
  
  emobj_pi_gr <- function(beta_pi){
    nsamp <- nrow(X)
    return(-t(X)%*%w/nsamp + t(X)%*%(1/(1+exp(-X%*%beta_pi)))/nsamp)
  }
  
  res <- optim(initl,emobj_pi,emobj_pi_gr,method = "BFGS")
  return(res$par)
}

normmixEM_contvar <- function(x_pi,x_m1,x_m2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_sgm21,initl_sgm22,maxiter=10000,epsilon=0.0001){ 
  beta_pi_old <- initl_beta_pi 
  beta_m1_old <- initl_beta_m1 
  beta_m2_old <- initl_beta_m2 
  sgm21_old <- initl_sgm21 
  sgm22_old <- initl_sgm22 
  w_base <- rep(0,length(y))
  l <- 0
  w_diff <- epsilon 
  
  while(l < maxiter & w_diff >= epsilon){
    l <- l + 1
    em_pi <- 1/(1+exp(-c(x_pi%*%beta_pi_old)))
    if(length(beta_m1_old)==1){em_m1 <- beta_m1_old}
    else{em_m1 <- c(x_m1%*%beta_m1_old)}
    if(length(beta_m2_old)==1){em_m2 <- beta_m2_old}
    else{em_m2 <- c(x_m2%*%beta_m2_old)}
    em_sgm21 <- sgm21_old
    em_sgm22 <- sgm22_old 
    w_update <- (em_pi*dnorm(y,mean=em_m1,sd=sqrt(em_sgm21)))/(em_pi*dnorm(y,mean=em_m1,sd=sqrt(em_sgm21))+(1-em_pi)*dnorm(y,mean=em_m2,sd=sqrt(em_sgm22)))
    beta_pi_new <- emopt_pi(beta_pi_old,x_pi,w_update)
    if(length(beta_m1_old)==1){beta_m1_new <- c(solve(t(x_m1)%*%diag(w_update)%*%x_m1)%*%t(x_m1)%*%diag(w_update)%*%y)
    sgm21_new <- sum(w_update*(y-beta_m1_new)^2)/sum(w_update)}
    else{beta_m1_new <- c(solve(t(x_m1)%*%diag(w_update)%*%x_m1)%*%t(x_m1)%*%diag(w_update)%*%y)
    sgm21_new <- sum(w_update*(y-x_m1%*%beta_m1_new)^2)/sum(w_update)}
    if(length(beta_m2_old)==1){beta_m2_new <- c(solve(t(x_m2)%*%diag(1-w_update)%*%x_m2)%*%t(x_m2)%*%diag(1-w_update)%*%y)
    sgm22_new <- sum((1-w_update)*(y-beta_m2_new)^2)/sum((1-w_update))}
    else{beta_m2_new <- c(solve(t(x_m2)%*%diag(1-w_update)%*%x_m2)%*%t(x_m2)%*%diag(1-w_update)%*%y)
    sgm22_new <- sum((1-w_update)*(y-x_m2%*%beta_m2_new)^2)/sum((1-w_update))}
    
    beta_pi_old <- beta_pi_new
    beta_m1_old <- beta_m1_new
    beta_m2_old <- beta_m2_new
    sgm21_old <- sgm21_new
    sgm22_old <- sgm22_new
    w_diff <- max(abs(w_update-w_base))
    w_base <- w_update
    
    #print(l)
    #print(beta_m1_new)
    
  }
  if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(beta_pi=beta_pi_old,beta_m1=beta_m1_old,beta_m2=beta_m2_old,sgm21=sgm21_old,sgm22=sgm22_old,iter=l))
}

loglik_mix2 <- function(y,lik.pi,lik.m1,lik.m2,lik.sgm21,lik.sgm22){
  return(sum(log(lik.pi*dnorm(y,mean=lik.m1,sd=sqrt(lik.sgm21))+(1-lik.pi)*dnorm(y,mean=lik.m2,sd=sqrt(lik.sgm21)))))
}

rnpnormmix <- function(p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}

pi_1 <- function(x){
  return(1/(1+exp(0.8-2*x))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(3 -0*x )   ###### beta = 0 
}

m_2 <- function(x){
  return(6 +2.5*x )
}

sigma_1 <- function(x){
  return(0.3*exp(1.5*x))
}

sigma_2 <- function(x){
  return(0.6*exp(-0.5*x))
}

n <- 200 
n_sim <- 500
n_bs <- 500
set.seed(599)
x <- runif(n)
p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(rep(0.7,n),rep(0.4,n))

X_pi <- X_m1_h1 <- X_m2 <- cbind(rep(1,n),x)
X_m1_h0 <- rep(1,n)

initl_beta_pi <- c(-0.8,2)
initl_beta_m1_h0 <- 3
initl_beta_m1_h1 <- c(3,0)
initl_beta_m2 <- c(6,2.5)
initl_sgm21 <- 0.49
initl_sgm22 <- 0.16
#### define a function for parallel simulation 
sim_mixreg_hp_mean <- function(j){
  library(mixtools)
  y <- rnpnormmix(p,mu,sgm)
  
  outsp_h1 <- normmixEM_contvar(X_pi,X_m1_h1,X_m2,y,initl_beta_pi,initl_beta_m1_h1,initl_beta_m2,initl_sgm21,initl_sgm22)
  outsp_h0 <- normmixEM_contvar(X_pi,X_m1_h0,X_m2,y,initl_beta_pi,initl_beta_m1_h0,initl_beta_m2,initl_sgm21,initl_sgm22)
  
  lik.pi.h1 <- 1/(1+exp(-c(X_pi%*%outsp_h1$beta_pi)))
  lik.pi.h0 <- 1/(1+exp(-c(X_pi%*%outsp_h0$beta_pi)))
  lik.m1.h1 <- c(X_m1_h1%*%outsp_h1$beta_m1)
  lik.m2.h1 <- c(X_m2%*%outsp_h1$beta_m2)
  lik.m1.h0 <- rep(outsp_h0$beta_m1,n)
  lik.m2.h0 <- c(X_m2%*%outsp_h0$beta_m2)
  T_lr <- 2*(loglik_mix2(y,lik.pi.h1,lik.m1.h1,lik.m2.h1,outsp_h1$sgm21,outsp_h1$sgm22)-loglik_mix2(y,lik.pi.h0,lik.m1.h0,lik.m2.h0,outsp_h0$sgm21,outsp_h0$sgm22))
  T_lrs <- rep(0,n_bs)
  
  p.mle <- cbind(lik.pi.h0,1-lik.pi.h0)
  mu.mle <- cbind(lik.m1.h0,lik.m2.h0)
  sgm.mle <- sqrt(cbind(rep(outsp_h0$sgm21,n),rep(outsp_h0$sgm22,n)))
  
  for(i_bs in 1:n_bs){
    y.bs <- rnpnormmix(p.mle,mu.mle,sgm.mle)
    
    outsp_h1_bs <- normmixEM_contvar(X_pi,X_m1_h1,X_m2,y.bs,initl_beta_pi,initl_beta_m1_h1,initl_beta_m2,initl_sgm21,initl_sgm22)
    outsp_h0_bs <- normmixEM_contvar(X_pi,X_m1_h0,X_m2,y.bs,initl_beta_pi,initl_beta_m1_h0,initl_beta_m2,initl_sgm21,initl_sgm22)
    
    lik.pi.h1_bs <- 1/(1+exp(-c(X_pi%*%outsp_h1_bs$beta_pi)))
    lik.pi.h0_bs <- 1/(1+exp(-c(X_pi%*%outsp_h0_bs$beta_pi)))
    lik.m1.h1_bs <- c(X_m1_h1%*%outsp_h1_bs$beta_m1)
    lik.m2.h1_bs <- c(X_m2%*%outsp_h1_bs$beta_m2)
    lik.m1.h0_bs <- rep(outsp_h0_bs$beta_m1,n)
    lik.m2.h0_bs <- c(X_m2%*%outsp_h0_bs$beta_m2)
    T_lrs[i_bs] <- 2*(loglik_mix2(y.bs,lik.pi.h1_bs,lik.m1.h1_bs,lik.m2.h1_bs,outsp_h1_bs$sgm21,outsp_h1_bs$sgm22)-
                        loglik_mix2(y.bs,lik.pi.h0_bs,lik.m1.h0_bs,lik.m2.h0_bs,outsp_h0_bs$sgm21,outsp_h0_bs$sgm22))
    
  }
  if(T_lr>quantile(T_lrs,0.95)) {rej <- 1}   #### alpha = 0.05
  else{rej <- 0}
  return(list(reject=rej))
}

the_cluster <- makeCluster(16)
clusterSetRNGStream(the_cluster,27)
clusterExport(the_cluster,c("emopt_pi","normmixEM_contvar","loglik_mix2","rnpnormmix","n","n_bs","p","mu","sgm","X_pi","X_m1_h1","X_m2","X_m1_h0","initl_beta_pi",
                            "initl_beta_m1_h0","initl_beta_m1_h1","initl_beta_m2","initl_sgm21","initl_sgm22"))
sim_hp <- parLapply(the_cluster,1:n_sim,sim_mixreg_hp_mean)
stopCluster(the_cluster)

reject <- sum(unlist(sim_hp))/n_sim
print(reject)



