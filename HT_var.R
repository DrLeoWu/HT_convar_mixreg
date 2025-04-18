##################################################################################################################
##################################################################################################################
####################   Simulation investigating the hypothesis testing on variances of mixture of regressions 
##################################################################################################################
##################################################################################################################
#### load libraries
library(mixtools)
library(sn)
library(matrixcalc)
library(Matrix)
library(ks)
library(pracma)
library(interp)
library(parallel)
########################################################################################################################
#### define functions 
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

#### define a function to get high dimensional ANOVA statistics
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

####  define a function to perform pseudoclass for two component mixture regression
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
  return(2 -2*x )   ###### beta = 0 
}

m_2 <- function(x){
  return(6 +3*x )
}

#sigma_1 <- function(x){
#  return(0.1+1*exp(-1.2+1.5*x))
#}
#sigma_1 <- function(x){
#  return(0.1+1*(x-0.5)^2)
#}

#sigma_1 <- function(x){
#  return(0.1+1*abs(sin(2*pi*x)))
#}
sigma_1 <- function(x){
  return(0.6+0.5*(x<0.7&x>0.3))
}
#sigma_1 <- 0.7
sigma_2 <- 0.4

n <- 200 
n_sim <- 500
n_bs <- 500
set.seed(599)
x <- runif(n)
p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),rep(sigma_2,n))

X_pi <- X_m1 <- X_m2 <- cbind(rep(1,n),x)

initl_beta_pi <- c(-0.8,2)
initl_beta_m1 <- c(3,-2)
initl_beta_m2 <- c(6,2.5)
initl_sgm21 <- 0.49
initl_sgm22 <- 0.16
#### define a function for parallel simulation 
sim_mixreg_hp_var <- function(j){
  library(mixtools)
  y <- rnpnormmix(p,mu,sgm)
  out <- normmixEM_contvar(X_pi,X_m1,X_m2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_sgm21,initl_sgm22)
  
  pi_est <- 1/(1+exp(-c(X_pi%*%out$beta_pi)))
  m1_est <- c(X_m1%*%out$beta_m1)
  m2_est <- c(X_m2%*%out$beta_m2)
  sd1_est <- sqrt(out$sgm21)
  sd2_est <- sqrt(out$sgm22)
  ### find the posterior probability
  post_pi <- (pi_est*dnorm(y,mean=m1_est,sd=sd1_est))/(pi_est*dnorm(y,mean=m1_est,sd=sd1_est)+(1-pi_est)*dnorm(y,mean=m2_est,sd=sd2_est))
  
  PseudoClass_Draw <- PseudoClass(x,y,post_pi,out$beta_m1,out$beta_m2)
  ANOVA_Sta <- Aug_ANOVA_T(PseudoClass_Draw$x_class1,abs(PseudoClass_Draw$adj_resd1),m=9)
  ANOVA_Sta_bs <- rep(0,n_bs)
  
  p.mle <- cbind(pi_est,1-pi_est)
  mu.mle <- cbind(m1_est,m2_est)
  sgm.mle <- cbind(rep(sd1_est,n),rep(sd2_est,n))
  ### start bootstrap
  for(i_bs in 1:n_bs){
    
    y.bs <- rnpnormmix(p.mle,mu.mle,sgm.mle)
    out_bs <- normmixEM_contvar(X_pi,X_m1,X_m2,y.bs,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_sgm21,initl_sgm22)
    
    pi_est_bs <- 1/(1+exp(-c(X_pi%*%out_bs$beta_pi)))
    m1_est_bs <- c(X_m1%*%out_bs$beta_m1)
    m2_est_bs <- c(X_m2%*%out_bs$beta_m2)
    sd1_est_bs <- sqrt(out_bs$sgm21)
    sd2_est_bs <- sqrt(out_bs$sgm22)
    ### find the posterior probability
    post_pi_bs <- (pi_est_bs*dnorm(y.bs,mean=m1_est_bs,sd=sd1_est_bs))/(pi_est_bs*dnorm(y.bs,mean=m1_est_bs,sd=sd1_est_bs)+(1-pi_est_bs)*dnorm(y.bs,mean=m2_est_bs,sd=sd2_est_bs))
    
    PseudoClass_Draw_bs <- PseudoClass(x,y.bs,post_pi_bs,out_bs$beta_m1,out_bs$beta_m2)
    ANOVA_Sta_bs[i_bs] <- Aug_ANOVA_T(PseudoClass_Draw_bs$x_class1,abs(PseudoClass_Draw_bs$adj_resd1),m=9)
    
  }
  
  if(ANOVA_Sta>quantile(ANOVA_Sta_bs,0.95)) {rej <- 1}   #### alpha = 0.05
  else {rej <- 0}
  return(list(reject=rej))
}

the_cluster <- makeCluster(16)
clusterSetRNGStream(the_cluster,36)
clusterExport(the_cluster,c("emopt_pi","normmixEM_contvar","Aug_ANOVA_T","PseudoClass","rnpnormmix","n","n_bs","p","mu","sgm","x","X_pi","X_m1","X_m2","initl_beta_pi",
                            "initl_beta_m1","initl_beta_m2","initl_sgm21","initl_sgm22"))
sim_hp <- parLapply(the_cluster,1:n_sim,sim_mixreg_hp_var)
stopCluster(the_cluster)

reject <- sum(unlist(sim_hp))/n_sim
print(reject)

