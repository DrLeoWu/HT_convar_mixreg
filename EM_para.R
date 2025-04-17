
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


emopt_wnorm <- function(initl_beta_m,initl_beta_sig2,X_m,X_sig2,y,w,maxiter=10000,epsilon=0.0001){  
  # this function solve the optimization problem of weighted log normal likelihood
  nsamp <- nrow(X_m)
  beta1_old <- initl_beta_m
  beta2_old <- initl_beta_sig2
  #w_new <- diag(w*exp(-X_sig2%*%beta2))
  l <- 0
  beta_diff <- epsilon 
  while(l < maxiter & beta_diff >= epsilon){
    l <- l + 1
    w_new <- diag(c(w*exp(-X_sig2%*%beta2_old)))
    beta1_new <- c(solve(t(X_m)%*%w_new%*%X_m)%*%t(X_m)%*%w_new%*%y)
    
    emobj_beta_sig2 <- function(beta_sig2){
      return(sum(w*(y-X_m%*%beta1_new)^2*exp(-X_sig2%*%beta_sig2)) + c(t(w)%*%(X_sig2%*%beta_sig2)))
    }
    
    emobj_beta_sig2_gr <- function(beta_sig2){
      return(-t(X_sig2)%*%(w*(y-X_m%*%beta1_new)^2*exp(-X_sig2%*%beta_sig2)) + t(X_sig2)%*%w)
    }
    
    rep2 <- optim(initl_beta_sig2,emobj_beta_sig2,emobj_beta_sig2_gr,method = "BFGS")
    beta2_new <- rep2$par
    beta_diff <- max(abs(c(beta1_new-beta1_old,beta2_new-beta2_old)))
    beta1_old <- beta1_new
    beta2_old <- beta2_new
  }
  
  return(list(beta_m=beta1_old,beta_sig2=beta2_old))
}




normmixEM <- function(x_pi,x_m,x_sgm2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22,maxiter=10000,epsilon=0.0001){ #x is n by p, and initl_beta is p by 5 (only for two clusters)
  beta_pi_old <- initl_beta_pi 
  beta_m1_old <- initl_beta_m1 
  beta_m2_old <- initl_beta_m2 
  beta_sgm21_old <- initl_beta_sgm21 
  beta_sgm22_old <- initl_beta_sgm22 
  w_base <- rep(0,length(y))
  l <- 0
  w_diff <- epsilon 
  
  while(l < maxiter & w_diff >= epsilon){
    l <- l + 1
    em_pi <- 1/(1+exp(-c(x_pi%*%beta_pi_old)))
    em_m1 <- c(x_m%*%beta_m1_old)
    em_m2 <- c(x_m%*%beta_m2_old)
    em_sgm21 <- exp(c(x_sgm2%*%beta_sgm21_old))
    em_sgm22 <- exp(c(x_sgm2%*%beta_sgm22_old))
    w_update <- (em_pi*dnorm(y,mean=em_m1,sd=sqrt(em_sgm21)))/(em_pi*dnorm(y,mean=em_m1,sd=sqrt(em_sgm21))+(1-em_pi)*dnorm(y,mean=em_m2,sd=sqrt(em_sgm22)))
    beta_pi_new <- emopt_pi(beta_pi_old,x_pi,w_update)
    normal1 <- emopt_wnorm(beta_m1_old,beta_sgm21_old,x_m,x_sgm2,y,w_update)
    normal2 <- emopt_wnorm(beta_m2_old,beta_sgm22_old,x_m,x_sgm2,y,1-w_update)
    beta_m1_new <- normal1$beta_m
    beta_m2_new <- normal2$beta_m
    beta_sgm21_new <- normal1$beta_sig2
    beta_sgm22_new <- normal2$beta_sig2
    
    beta_pi_old <- beta_pi_new
    beta_m1_old <- beta_m1_new
    beta_m2_old <- beta_m2_new
    beta_sgm21_old <- beta_sgm21_new
    beta_sgm22_old <- beta_sgm22_new
    w_diff <- max(abs(w_update-w_base))
    w_base <- w_update
    
    #print(l)
    #print(w_base)

  }
  if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(beta_pi=beta_pi_old,beta_m1=beta_m1_old,beta_m2=beta_m2_old,beta_sgm21=beta_sgm21_old,beta_sgm22=beta_sgm22_old,iter=l))
}

#X_reform <- cbind(rep(1,n),x,x^2,x[,1]*x[,2])
pi_1 <- function(x){
  return(exp(0.5*x[,1])/(1+exp(0.5*x[,2]))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  #return(6-sin(2*pi*x[,1])-sin(2*pi*x[,2]))
  return(6-x[,1]-2*x[,2])
}

m_2 <- function(x){
  #return(cos(1*pi*x[,1]) + cos(4*pi*x[,2]))
  return(-1+2*x[,1]+x[,2])
}

sigma_1 <- function(x){
  return(0.6*exp(0.1*x[,1]+0.4*x[,2]))
}

sigma_2 <- function(x){
  return(0.5*exp(-0.1*x[,1]-0.2*x[,2]))
}

rnpnormmix <- function(x,p,m,sgm){
  y <- NULL
  for (i in 1:nrow(x)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}

n <- 100 
set.seed(599)
#set.seed(2468)
x <- cbind(runif(n),runif(n))
#x <- cbind(runif(n),rnorm(n,2))
#p <- cbind(pi_1(x),pi_2(x))
p <- cbind(rep(0.3,n),rep(0.7,n))
mu <- cbind(m_1(x),m_2(x))
#sgm <- cbind(sigma_1(x),sigma_2(x))
sgm <- cbind(rep(0.7,n),rep(0.4,n))
y <- rnpnormmix(x,p,mu,sgm)

rmt <- regmix(x,y,nclust=2)
print(rmt$coef)
print(rmt$var)
print(rmt$eps)

X_pi <- X_sig2 <- cbind(rep(1,n),x)
#X_m <- cbind(rep(1,n),x)
X_m <- cbind(rep(1,n),x,x^2,x[,1]*x[,2])
#X_m <- cbind(rep(1,n1),x,x^2,x^3,x^4)
initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,3)
initl_beta_m1 <- initl_beta_m2 <- rep(0.5,6)
test <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)

#####################################################################################################################################
#####################################################################################################################################
################### The fowllowing function perform EM for 2 component mixture with constant variance 
normmixEM_contvar <- function(x_pi,x_m,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_sgm21,initl_sgm22,maxiter=10000,epsilon=0.0001){ #x is n by p, and initl_beta is p by 5 (only for two clusters)
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
    em_m1 <- c(x_m%*%beta_m1_old)
    em_m2 <- c(x_m%*%beta_m2_old)
    em_sgm21 <- sgm21_old
    em_sgm22 <- sgm22_old 
    w_update <- (em_pi*dnorm(y,mean=em_m1,sd=sqrt(em_sgm21)))/(em_pi*dnorm(y,mean=em_m1,sd=sqrt(em_sgm21))+(1-em_pi)*dnorm(y,mean=em_m2,sd=sqrt(em_sgm22)))
    beta_pi_new <- emopt_pi(beta_pi_old,x_pi,w_update)
    #normal1 <- emopt_wnorm(beta_m1_old,beta_sgm21_old,x_m,x_sgm2,y,w_update)
    #normal2 <- emopt_wnorm(beta_m2_old,beta_sgm22_old,x_m,x_sgm2,y,1-w_update)
    beta_m1_new <- c(solve(t(x_m)%*%diag(w_update)%*%x_m)%*%t(x_m)%*%diag(w_update)%*%y)
    beta_m2_new <- c(solve(t(x_m)%*%diag(1-w_update)%*%x_m)%*%t(x_m)%*%diag(1-w_update)%*%y)
    sgm21_new <- sum(w_update*(y-x_m%*%beta_m1_new)^2)/sum(w_update)
    sgm22_new <- sum((1-w_update)*(y-x_m%*%beta_m2_new)^2)/sum((1-w_update))
    
    beta_pi_old <- beta_pi_new
    beta_m1_old <- beta_m1_new
    beta_m2_old <- beta_m2_new
    sgm21_old <- sgm21_new
    sgm22_old <- sgm22_new
    w_diff <- max(abs(w_update-w_base))
    w_base <- w_update
    
    #print(l)
    #print(w_base)
    
  }
  if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(beta_pi=beta_pi_old,beta_m1=beta_m1_old,beta_m2=beta_m2_old,sgm21=sgm21_old,sgm22=sgm22_old,iter=l))
}

#####################################################################################################################################
################### The fowllowing function perform EM for 2 component mixture with constant variance and means
normmixEM_cont_mv <- function(x_pi,y,initl_beta_pi,initl_m1,initl_m2,initl_sgm21,initl_sgm22,maxiter=10000,epsilon=0.0001){ #x is n by p, and initl_beta is p by 5 (only for two clusters)
  beta_pi_old <- initl_beta_pi 
  m1_old <- initl_m1 
  m2_old <- initl_m2 
  sgm21_old <- initl_sgm21 
  sgm22_old <- initl_sgm22 
  w_base <- rep(0,length(y))
  l <- 0
  w_diff <- epsilon 
  
  while(l < maxiter & w_diff >= epsilon){
    l <- l + 1
    em_pi <- 1/(1+exp(-c(x_pi%*%beta_pi_old)))
    em_m1 <- m1_old
    em_m2 <- m2_old
    em_sgm21 <- sgm21_old
    em_sgm22 <- sgm22_old 
    w_update <- (em_pi*dnorm(y,mean=em_m1,sd=sqrt(em_sgm21)))/(em_pi*dnorm(y,mean=em_m1,sd=sqrt(em_sgm21))+(1-em_pi)*dnorm(y,mean=em_m2,sd=sqrt(em_sgm22)))
    beta_pi_new <- emopt_pi(beta_pi_old,x_pi,w_update)
    #normal1 <- emopt_wnorm(beta_m1_old,beta_sgm21_old,x_m,x_sgm2,y,w_update)
    #normal2 <- emopt_wnorm(beta_m2_old,beta_sgm22_old,x_m,x_sgm2,y,1-w_update)
    m1_new <- sum(w_update*y)/sum(w_update)
    m2_new <- sum((1-w_update)*y)/sum(1-w_update)
    sgm21_new <- sum(w_update*(y-m1_new)^2)/sum(w_update)
    sgm22_new <- sum((1-w_update)*(y-m2_new)^2)/sum((1-w_update))
    
    beta_pi_old <- beta_pi_new
    m1_old <- m1_new
    m2_old <- m2_new
    sgm21_old <- sgm21_new
    sgm22_old <- sgm22_new
    w_diff <- max(abs(w_update-w_base))
    w_base <- w_update
    
    #print(l)
    #print(w_base)
    
  }
  if (l >= maxiter) {print("The algorithm does not converge")}
  return(list(beta_pi=beta_pi_old,m1=m1_old,m2=m2_old,sgm21=sgm21_old,sgm22=sgm22_old,iter=l))
}

#####################################################################################################################################
#####################################################################################################################################
################### Hypothesis testing for the normal mixture regression 

loglik_mix2 <- function(y,lik.pi,lik.m1,lik.m2,lik.sgm21,lik.sgm22){
  #lik.pi <- 1/(1+exp(-c(xd%*%beta_pi)))
  #lik.m1 <- c(xd%*%beta_m1)
  #lik.m2 <- c(xd%*%beta_m2)
  #lik.sgm21 <- exp(c(xd%*%beta_sgm21))
  #lik.sgm22 <- exp(c(xd%*%beta_sgm22))

  return(sum(log(lik.pi*dnorm(y,mean=lik.m1,sd=sqrt(lik.sgm21))+(1-lik.pi)*dnorm(y,mean=lik.m2,sd=sqrt(lik.sgm21)))))
}

pi_1 <- function(x){
  return(exp(2-4*x)/(1+exp(2-4*x))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  #return(6-sin(2*pi*x[,1])-sin(2*pi*x[,2]))
  return(5-x)
}

m_2 <- function(x){
  #return(cos(1*pi*x[,1]) + cos(4*pi*x[,2]))
  return(1+2*x)
}



sigma_1 <- function(x){
  return(exp(-0.5+0.5*x))
}

sigma_2 <- function(x){
  return(exp(-0.1*x[,1]-0.2*x[,2]))
}

rnpnormmix <- function(p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}


n <- 200 
n_sim <- 500
set.seed(599)
x <- runif(n)
p <- cbind(pi_1(x),pi_2(x))
#mu <- cbind(m_1(x),m_2(x))
mu <- cbind(rep(4,n),rep(1,n))
sgm <- cbind(rep(0.7,n),rep(0.4,n))

#rmt <- regmix(x,y,nclust=2)
#print(rmt$coef)
#print(rmt$var)
#print(rmt$eps)

X_pi <- X_m <- cbind(rep(1,n),x)

initl_beta_pi <- c(1,-1)
initl_beta_m1 <- c(3,0.1)
initl_beta_m2 <- c(0.5,0.1)

rej <- 0
n_bs <- 500

for(i_sim in 1:n_sim){
  y <- rnpnormmix(p,mu,sgm)
  
  outsp_h1 <- normmixEM_contvar(X_pi,X_m,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_sgm21=0.5,initl_sgm22=0.5,maxiter=10000,epsilon=0.0001)
  outsp_h0 <- normmixEM_cont_mv(X_pi,y,initl_beta_pi,initl_m1=3,initl_m2=0.5,initl_sgm21=0.5,initl_sgm22=0.5,maxiter=10000,epsilon=0.0001)
  
  lik.pi.h1 <- 1/(1+exp(-c(X_pi%*%outsp_h1$beta_pi)))
  lik.pi.h0 <- 1/(1+exp(-c(X_pi%*%outsp_h0$beta_pi)))
  lik.m1.h1 <- c(X_m%*%outsp_h1$beta_m1)
  lik.m2.h1 <- c(X_m%*%outsp_h1$beta_m2)
  T_lr <- 2*(loglik_mix2(y,lik.pi.h1,lik.m1.h1,lik.m2.h1,outsp_h1$sgm21,outsp_h1$sgm22)-loglik_mix2(y,lik.pi.h0,outsp_h0$m1,outsp_h0$m2,outsp_h0$sgm21,outsp_h0$sgm22))
  T_lrs <- rep(0,n_bs)
  for(i_bs in 1:n_bs){
    
    p.mle <- cbind(lik.pi.h0,1-lik.pi.h0)
    mu.mle <- cbind(rep(outsp_h0$m1,n),rep(outsp_h0$m2,n))
    #sgm <- cbind(sigma_1(x),sigma_2(x))
    sgm.mle <- sqrt(cbind(rep(outsp_h0$sgm21,n),rep(outsp_h0$sgm22,n)))
    y.bs <- rnpnormmix(p.mle,mu.mle,sgm.mle)
    
    outsp_h1_bs <- normmixEM_contvar(X_pi,X_m,y.bs,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_sgm21=0.5,initl_sgm22=0.5,maxiter=10000,epsilon=0.0001)
    outsp_h0_bs <- normmixEM_cont_mv(X_pi,y.bs,initl_beta_pi,initl_m1=3,initl_m2=0.5,initl_sgm21=0.5,initl_sgm22=0.5,maxiter=10000,epsilon=0.0001)
    
    lik.pi.h1_bs <- 1/(1+exp(-c(X_pi%*%outsp_h1_bs$beta_pi)))
    lik.pi.h0_bs <- 1/(1+exp(-c(X_pi%*%outsp_h0_bs$beta_pi)))
    lik.m1.h1_bs <- c(X_m%*%outsp_h1_bs$beta_m1)
    lik.m2.h1_bs <- c(X_m%*%outsp_h1_bs$beta_m2)
    T_lrs[i_bs] <- 2*(loglik_mix2(y.bs,lik.pi.h1_bs,lik.m1.h1_bs,lik.m2.h1_bs,outsp_h1_bs$sgm21,outsp_h1_bs$sgm22)
                      -loglik_mix2(y.bs,lik.pi.h0_bs,outsp_h0_bs$m1,outsp_h0_bs$m2,outsp_h0_bs$sgm21,outsp_h0_bs$sgm22))
    
  }
  if(T_lr>quantile(T_lrs,0.95)) {rej <- rej + 1}
  print(i_sim)
}



#plot(density(T_lrs))


#test2 <- hmeEM(y,x,beta=cbind(initl_beta_m1,initl_beta_m2))
#set.seed(100)
#em.out <- regmixEM(Equivalence, NO)
#hme.out <- hmeEM(Equivalence, NO, beta = em.out$beta)
pi_1 <- function(x){
  return(exp(2-4*x)/(1+exp(2-4*x))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  #return(6-sin(2*pi*x[,1])-sin(2*pi*x[,2]))
  return(5-x)
}

m_2 <- function(x){
  #return(cos(1*pi*x[,1]) + cos(4*pi*x[,2]))
  return(1+2*x)
}



sigma_1 <- function(x){
  return(exp(-0.5+0.6*x))
}

sigma_2 <- function(x){
  return(exp(-0.6-0.3*x))
}

n <- 200 
n_sim <- 500
set.seed(599)
x <- runif(n)
p <- cbind(pi_1(x),pi_2(x))
#mu <- cbind(m_1(x),m_2(x))
mu <- cbind(rep(4,n),rep(1,n))
sgm <- cbind(sigma_1(x),sigma_2(x))

#rmt <- regmix(x,y,nclust=2)
#print(rmt$coef)
#print(rmt$var)
#print(rmt$eps)

X_pi <- X_m <- cbind(rep(1,n),x)

initl_beta_pi <- c(1,-1)
initl_beta_m1 <- c(3,0.1)
initl_beta_m2 <- c(0.5,0.1)

rej <- 0
n_bs <- 500

for(i_sim in 1:n_sim){
  y <- rnpnormmix(p,mu,sgm)
  
  outsp_h1 <- normmixEM_contvar(X_pi,X_m,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_sgm21=0.5,initl_sgm22=0.5,maxiter=10000,epsilon=0.0001)
  outsp_h0 <- normmixEM_cont_mv(X_pi,y,initl_beta_pi,initl_m1=3,initl_m2=0.5,initl_sgm21=0.5,initl_sgm22=0.5,maxiter=10000,epsilon=0.0001)
  
  lik.pi.h1 <- 1/(1+exp(-c(X_pi%*%outsp_h1$beta_pi)))
  lik.pi.h0 <- 1/(1+exp(-c(X_pi%*%outsp_h0$beta_pi)))
  lik.m1.h1 <- c(X_m%*%outsp_h1$beta_m1)
  lik.m2.h1 <- c(X_m%*%outsp_h1$beta_m2)
  T_lr <- 2*(loglik_mix2(y,lik.pi.h1,lik.m1.h1,lik.m2.h1,outsp_h1$sgm21,outsp_h1$sgm22)-loglik_mix2(y,lik.pi.h0,outsp_h0$m1,outsp_h0$m2,outsp_h0$sgm21,outsp_h0$sgm22))
  T_lrs <- rep(0,n_bs)
  for(i_bs in 1:n_bs){
    
    p.mle <- cbind(lik.pi.h0,1-lik.pi.h0)
    mu.mle <- cbind(rep(outsp_h0$m1,n),rep(outsp_h0$m2,n))
    #sgm <- cbind(sigma_1(x),sigma_2(x))
    sgm.mle <- sqrt(cbind(rep(outsp_h0$sgm21,n),rep(outsp_h0$sgm22,n)))
    y.bs <- rnpnormmix(p.mle,mu.mle,sgm.mle)
    
    outsp_h1_bs <- normmixEM_contvar(X_pi,X_m,y.bs,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_sgm21=0.5,initl_sgm22=0.5,maxiter=10000,epsilon=0.0001)
    outsp_h0_bs <- normmixEM_cont_mv(X_pi,y.bs,initl_beta_pi,initl_m1=3,initl_m2=0.5,initl_sgm21=0.5,initl_sgm22=0.5,maxiter=10000,epsilon=0.0001)
    
    lik.pi.h1_bs <- 1/(1+exp(-c(X_pi%*%outsp_h1_bs$beta_pi)))
    lik.pi.h0_bs <- 1/(1+exp(-c(X_pi%*%outsp_h0_bs$beta_pi)))
    lik.m1.h1_bs <- c(X_m%*%outsp_h1_bs$beta_m1)
    lik.m2.h1_bs <- c(X_m%*%outsp_h1_bs$beta_m2)
    T_lrs[i_bs] <- 2*(loglik_mix2(y.bs,lik.pi.h1_bs,lik.m1.h1_bs,lik.m2.h1_bs,outsp_h1_bs$sgm21,outsp_h1_bs$sgm22)
                      -loglik_mix2(y.bs,lik.pi.h0_bs,outsp_h0_bs$m1,outsp_h0_bs$m2,outsp_h0_bs$sgm21,outsp_h0_bs$sgm22))
    
  }
  if(T_lr>quantile(T_lrs,0.95)) {rej <- rej + 1}
  print(i_sim)
}






