#########################################################################################################################################
########  Simulation studies to evaluate the performance of the proposed hypothesis test on variance constancy for two-component mixtures
#########################################################################################################################################


#### Load libraries
library(mixtools)
library(sn)
library(matrixcalc)
library(Matrix)
library(ks)
library(pracma)
library(interp)
library(parallel)

########################################################################################################################
#################### Define a function to get high dimensional ANOVA statistics
Aug_ANOVA_T <- function(x,y,m){#  x and y are univariate and m (an odd number) is the window size
  n_cell <- length(x)
  y_order <- y[order(x)]
  MSE <- 0
  sum_cell <- mean_cell <- rep(0,n_cell)
  for(i in 1:n_cell){
    
    if (i < (m+1)/2 ){
      mean_cell[i] <- mean(y_order[1:m])
      sum_cell[i] <- sum(y_order[1:m])
      MSE <- MSE + sum((y_order[1:m]-mean_cell[i])^2)
    } else if(i > (n_cell-(m-1)/2)){
      mean_cell[i] <- mean(y_order[(n_cell-m+1):n_cell])
      sum_cell[i] <- sum(y_order[(n_cell-m+1):n_cell])
      MSE <- MSE + sum((y_order[(n_cell-m+1):n_cell]-mean_cell[i])^2)
    } else {
      mean_cell[i] <- mean(y_order[(i-(m-1)/2):(i+(m-1)/2)])
      sum_cell[i] <- sum(y_order[(i-(m-1)/2):(i+(m-1)/2)])
      MSE <- MSE + sum((y_order[(i-(m-1)/2):(i+(m-1)/2)]-mean_cell[i])^2)
    }
    
  }
  MSE <- MSE/(n_cell*(m-1))
  MST <- m*sum((mean_cell - sum(sum_cell)/(n_cell*m))^2)/(n_cell-1)
  StaT <- sqrt(n_cell)*(MST-MSE)
  return(StaT)
}
################  define a function to perform pseudoclass for two component mixture regression
PseudoClass <- function(x,y,p1,beta1,beta2){# x and y are univariate, and p1 is the posterior probability
  randraw <- rbinom(length(y),1,p1)
  y_c1 <- y[randraw==1]
  y_c2 <- y[randraw==0]
  x_c1 <- x[randraw==1]
  x_c2 <- x[randraw==0]
  xd1 <- cbind(rep(1,length(x_c1)),x_c1)
  xd2 <- cbind(rep(1,length(x_c2)),x_c2)
  adj_resd1 <- y_c1 - xd1%*%beta1
  adj_resd2 <- y_c2 - xd2%*%beta2
  return(list(y_class1=y_c1,y_class2=y_c2,x_class1=x_c1,x_class2=x_c2,adj_resd1=adj_resd1,adj_resd2=adj_resd2))
}


################  define functions to perform EM algorithm for two component mixture regression with covariate dependent mixing proportions and constant variance
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

########################################################################################################################
#################### start simulation studies 

## define a function to generate samples from mixtures 
rnpnormmix <- function(p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}

######################################### setting 1 #######################################################################################

pi_1 <- function(x){
  return(exp(2-4*x)/(1+exp(2-4*x))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(5-x)
}

m_2 <- function(x){
  return(1+2*x)
}

sigma_1 <- 0.7; sigma_2 <- 0.4

n <- 1000 
n_sim <- 1000
n_impt <- 1 # n_impt is the number of independent imputation
set.seed(489)
x <- runif(n)
p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(rep(sigma_1,n),rep(sigma_2,n))

true_beta_pi <- c(2,-4)
true_beta_m1 <- c(5,-1)
true_beta_m2 <- c(1,2)
x_m <- cbind(rep(1,n),x)
m1 <- c(x_m%*%true_beta_m1)
m2 <- c(x_m%*%true_beta_m2)
ANOVA_Sta <- rep(0,n_sim)

for(i in 1:n_sim){
  
  y <- rnpnormmix(p,mu,sgm)
  post_pi <- (pi_1(x)*dnorm(y,mean=m1,sd=sigma_1)/(pi_1(x)*dnorm(y,mean=m1,sd=sigma_1)+pi_2(x)*dnorm(y,mean=m2,sd=sigma_2)))
  ANOVA_Stas_impt <- rep(0,n_impt)
  for(k in 1:n_impt){  
    PseudoClass_Draw <- PseudoClass(x,y,post_pi,true_beta_m1,true_beta_m2)
    ANOVA_Stas_impt[k] <- Aug_ANOVA_T(PseudoClass_Draw$x_class1,abs(PseudoClass_Draw$adj_resd1),m=7)
  }
  ANOVA_Sta[i] <- mean(ANOVA_Stas_impt)
}

m <- 7
#Sta_sd <- var(abs(PseudoClass_Draw$adj_resd1))*sqrt(2*m*(2*m-1)/(3*(m-1)))
theo_sd <- sigma_1^2*(1-2/pi)*sqrt(2*m*(2*m-1)/(3*(m-1)))
plot(density(ANOVA_Sta),ylim=c(0,1))
lines(seq(-3,3,length=100),dnorm(seq(-3,3,length=100),mean=0,sd=theo_sd),type="l",col="red")
up_bound <- qnorm(0.975,mean=0,sd=theo_sd)
low_bound <- qnorm(0.025,mean=0,sd=theo_sd)
empir_alpha <- length(ANOVA_Sta[ANOVA_Sta>=up_bound | ANOVA_Sta<=low_bound])/n_sim
print(empir_alpha)

#plot(density(rnorm(5000)))
#lines(seq(-3,3,length=100),dnorm(seq(-3,3,length=100)),type="l",col="red")

######################################### setting 2 #######################################################################################

pi_1 <- function(x){
  return(exp(-1+2*x)/(1+exp(-1+2*x))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(3+2*x)
}

m_2 <- function(x){
  return(1-x)
}

sigma_1 <- 0.8; sigma_2 <- 0.5

n <- 1000 
n_sim <- 1000
n_impt <- 1 # n_impt is the number of independent imputation
set.seed(489)
x <- runif(n)
p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(rep(sigma_1,n),rep(sigma_2,n))

true_beta_pi <- c(-1,2)
true_beta_m1 <- c(3,2)
true_beta_m2 <- c(1,-1)
x_m <- cbind(rep(1,n),x)
m1 <- c(x_m%*%true_beta_m1)
m2 <- c(x_m%*%true_beta_m2)
ANOVA_Sta <- rep(0,n_sim)

for(i in 1:n_sim){
  
  y <- rnpnormmix(p,mu,sgm)
  post_pi <- (pi_1(x)*dnorm(y,mean=m1,sd=sigma_1)/(pi_1(x)*dnorm(y,mean=m1,sd=sigma_1)+pi_2(x)*dnorm(y,mean=m2,sd=sigma_2)))
  ANOVA_Stas_impt <- rep(0,n_impt)
  for(k in 1:n_impt){  
    PseudoClass_Draw <- PseudoClass(x,y,post_pi,true_beta_m1,true_beta_m2)
    ANOVA_Stas_impt[k] <- Aug_ANOVA_T(PseudoClass_Draw$x_class1,abs(PseudoClass_Draw$adj_resd1),m=7)
  }
  ANOVA_Sta[i] <- mean(ANOVA_Stas_impt)
}

m <- 7
#Sta_sd <- var(abs(PseudoClass_Draw$adj_resd1))*sqrt(2*m*(2*m-1)/(3*(m-1)))
theo_sd <- sigma_1^2*(1-2/pi)*sqrt(2*m*(2*m-1)/(3*(m-1)))
plot(density(ANOVA_Sta),ylim=c(0,1))
lines(seq(-3,3,length=100),dnorm(seq(-3,3,length=100),mean=0,sd=theo_sd),type="l",col="red")
up_bound <- qnorm(0.975,mean=0,sd=theo_sd)
low_bound <- qnorm(0.025,mean=0,sd=theo_sd)
empir_alpha <- length(ANOVA_Sta[ANOVA_Sta>=up_bound | ANOVA_Sta<=low_bound])/n_sim
print(empir_alpha)

######################################### setting 3 #######################################################################################
pi_1 <- function(x){
  return(exp(0+1*x)/(1+exp(0+1*x))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(3+1*x)
}

m_2 <- function(x){
  return(-2+3*x)
}

sigma_1 <- 0.6; sigma_2 <- 0.8

n <- 1000 
n_sim <- 1000
n_impt <- 1 # n_impt is the number of independent imputation
set.seed(489)
x <- runif(n)
p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(rep(sigma_1,n),rep(sigma_2,n))

true_beta_pi <- c(0,1)
true_beta_m1 <- c(3,1)
true_beta_m2 <- c(-2,3)
x_m <- cbind(rep(1,n),x)
m1 <- c(x_m%*%true_beta_m1)
m2 <- c(x_m%*%true_beta_m2)
ANOVA_Sta <- rep(0,n_sim)

for(i in 1:n_sim){
  
  y <- rnpnormmix(p,mu,sgm)
  post_pi <- (pi_1(x)*dnorm(y,mean=m1,sd=sigma_1)/(pi_1(x)*dnorm(y,mean=m1,sd=sigma_1)+pi_2(x)*dnorm(y,mean=m2,sd=sigma_2)))
  ANOVA_Stas_impt <- rep(0,n_impt)
  for(k in 1:n_impt){  
    PseudoClass_Draw <- PseudoClass(x,y,post_pi,true_beta_m1,true_beta_m2)
    ANOVA_Stas_impt[k] <- Aug_ANOVA_T(PseudoClass_Draw$x_class1,abs(PseudoClass_Draw$adj_resd1),m=7)
  }
  ANOVA_Sta[i] <- mean(ANOVA_Stas_impt)
}

m <- 7
#Sta_sd <- var(abs(PseudoClass_Draw$adj_resd1))*sqrt(2*m*(2*m-1)/(3*(m-1)))
theo_sd <- sigma_1^2*(1-2/pi)*sqrt(2*m*(2*m-1)/(3*(m-1)))
plot(density(ANOVA_Sta),ylim=c(0,1))
lines(seq(-3,3,length=100),dnorm(seq(-3,3,length=100),mean=0,sd=theo_sd),type="l",col="red")
up_bound <- qnorm(0.975,mean=0,sd=theo_sd)
low_bound <- qnorm(0.025,mean=0,sd=theo_sd)
empir_alpha <- length(ANOVA_Sta[ANOVA_Sta>=up_bound | ANOVA_Sta<=low_bound])/n_sim
print(empir_alpha)





