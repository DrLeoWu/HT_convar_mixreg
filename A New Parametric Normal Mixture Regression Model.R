
x <- seq(0,1,0.01)

pi_1 <- function(x){
  return(exp(0.5*x)/(1+exp(0.5*x)))
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(3-sin(2*pi*x))
}

m_2 <- function(x){
  return(cos(3*pi*x))
}

sigma_1 <- function(x){
  return(0.6*exp(0.5*x))
}

sigma_2 <- function(x){
  return(0.5*exp(-0.2*x))
}

plot(x,m_1(x),type="l")
lines(x,m_2(x),col=5)
plot(x,m_2(x),type="l",add=TRUE,col=3)

plot( x, m_1(x), type="l", col="red" )
par(new=TRUE)
plot( x, m_2(x), type="l", col="green" )



par( mfrow= c(3,2) )

plot(x,pis,type="l",ylab="pi(x)",ylim=c(0,1))
pis <- 1/(1+exp(0.4916863-1.1856913*x))
lines(x,pi_1(x),col="red")

y1 <- m_1(x)
y2 <- 4.618536 -39.185857*x +  75.716174*x^2 -17.960661*x^3 -20.846463*x^4

plot(x,y2,type="l",ylab="m1(x)")
lines(x,y1,col="red")


m2 <- 1.541159  -5.885048*x + 50.691695*x^2 -94.000114*x^3 +  46.595915*x^4
plot(x,m2,type="l",ylab="m2(x)")
lines(x,m_2(x),col="red")

sig2_1 <- exp(-0.5690471 + 0.8161908*x)
sig2_2 <- exp(0.5435652 -2.6158204*x)

plot(x,sig2_1,type="l",ylab="sigam2_1(x)",ylim=c(0,2))
lines(x,sigma_1(x)^2,col="red")

plot(x,sig2_2,type="l",ylab="sigam2_2(x)",ylim=c(0,2))
lines(x,sigma_2(x)^2,col="red")
##############################################################################################################
###################### test the validity of the parametric EM algorithm 

###### set the true functionals
set.seed(12345)
n <- 200
x <- runif(n)

pi_1 <- function(x){
  return(1/(1+exp(0.8-2*x))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(2 -1*x +  6*x^2 -2.5*x^3 -3.5*x^4)
  }

m_2 <- function(x){
  return(6 -1*x -  4*x^2 +1.5*x^3 +2.5*x^4)
  }

sigma_1 <- function(x){
  return(0.6*exp(0.5*x))
}

sigma_2 <- function(x){
  return(0.5*exp(-0.2*x))
}

rnpnormmix <- function(x,p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}


p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))
y <- rnpnormmix(x,p,mu,sgm)

X_pi <- X_sig2 <- cbind(rep(1,n),x)
#X_m <- cbind(rep(1,n),x,x^2,x[,1]*x[,2])
X_m <- cbind(rep(1,n),x,x^2,x^3,x^4)
initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,2)
initl_beta_m1 <- initl_beta_m2 <- rep(0.5,5)
test <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)

x00 <- seq(0,1,0.01)

pi_hat <- function(x00){
  return(1/(1+exp(cbind(rep(1,length(x00)),x00)%*%c(test$beta_pi))))
}

m1_hat <- function(x00){
  return(cbind(rep(1,length(x00)),x00,x00^2,x00^3,x00^4)%*%c(test$beta_m1))
}


m2_hat <- function(x00){
  return(cbind(rep(1,length(x00)),x00,x00^2,x00^3,x00^4)%*%c(test$beta_m2))
}

sigma21_hat <- function(x00){
  return(exp(cbind(rep(1,length(x00)),x00)%*%c(test$beta_sgm21)))
}

sigma22_hat <- function(x00){
  return(exp(cbind(rep(1,length(x00)),x00)%*%c(test$beta_sgm22)))
}



par( mfrow= c(3,2) )

plot(x00,pi_hat(x00),type="l",ylab="pi(x)",ylim=c(0,1))
lines(x00,pi_1(x00),col="red")

plot(x00,m1_hat(x00),type="l",ylab="m1(x)")
lines(x00,m_2(x00),col="red")

plot(x00,m2_hat(x00),type="l",ylab="m2(x)")
lines(x00,m_1(x00),col="red")

plot(x00,sigma22_hat(x00),type="l",ylab="sigma22",ylim=c(0.3,1))
lines(x00,sigma_1(x00)^2,col="red")

plot(x00,sigma21_hat(x00),type="l",ylab="sigma21",ylim=c(0.1,0.6))
lines(x00,sigma_2(x00)^2,col="red")

############################################################################################

ptheta1_1 <- exp(-X_pi%*%test$beta_pi)/(1+exp(-X_pi%*%test$beta_pi))^2*test$beta_pi[2]
ptheta2_1 <- exp(X_sig2%*%test$beta_sgm21)*test$beta_sgm21[2]

ptheta3_1 <- exp(X_sig2%*%test$beta_sgm22)*test$beta_sgm22[2] 

ptheta4_1 <- test$beta_m1[2] + 2*test$beta_m1[4]*x[,1] + test$beta_m1[6]*x[,2]

ptheta5_1 <- test$beta_m2[2] + 2*test$beta_m2[4]*x[,1] + test$beta_m2[6]*x[,2]

ptheta_1 <- cbind(x[,1],ptheta1_1,ptheta2_1,ptheta3_1,ptheta4_1,ptheta5_1)

ptheta1_2 <- exp(-X_pi%*%test$beta_pi)/(1+exp(-X_pi%*%test$beta_pi))^2*test$beta_pi[3]

ptheta2_2 <- exp(X_sig2%*%test$beta_sgm21)*test$beta_sgm21[3]

ptheta3_2 <- exp(X_sig2%*%test$beta_sgm22)*test$beta_sgm22[3]

ptheta4_2 <- test$beta_m1[3] + 2*test$beta_m1[5]*x[,1] + test$beta_m1[6]*x[,1]

ptheta5_2 <- test$beta_m2[3] + 2*test$beta_m2[5]*x[,1] + test$beta_m2[6]*x[,1]

ptheta_2 <- cbind(x[,2],ptheta1_2,ptheta2_2,ptheta3_2,ptheta4_2,ptheta5_2)

pi_intm <- exp(-X_pi%*%test$beta_pi)

ptheta1_11 <- test$beta_pi[2]^2*pi_intm*(-1+pi_intm)/(1+pi_intm)^3

ptheta2_11 <- exp(X_sig2%*%test$beta_sgm21)*test$beta_sgm21[2]^2

ptheta3_11 <- exp(X_sig2%*%test$beta_sgm22)*test$beta_sgm22[2]^2

ptheta4_11 <- rep(2*test$beta_m1[4],nrow(x))

ptheta5_11 <- rep(2*test$beta_m2[4],nrow(x))

ptheta_11 <- cbind(x[,1],ptheta1_11,ptheta2_11,ptheta3_11,ptheta4_11,ptheta5_11)

ptheta1_22 <- test$beta_pi[3]^2*pi_intm*(-1+pi_intm)/(1+pi_intm)^3

ptheta2_22 <- exp(X_sig2%*%test$beta_sgm21)*test$beta_sgm21[3]^2

ptheta3_22 <- exp(X_sig2%*%test$beta_sgm22)*test$beta_sgm22[3]^2

ptheta4_22 <- rep(2*test$beta_m1[5],nrow(x))

ptheta5_22 <- rep(2*test$beta_m2[5],nrow(x))

ptheta_22 <- cbind(x[,2],ptheta1_22,ptheta2_22,ptheta3_22,ptheta4_22,ptheta5_22)

ptheta1_12 <- test$beta_pi[2]*test$beta_pi[3]*pi_intm*(-1+pi_intm)/(1+pi_intm)^3

ptheta2_12 <- exp(X_sig2%*%test$beta_sgm21)*test$beta_sgm21[2]*test$beta_sgm21[3]

ptheta3_12 <- exp(X_sig2%*%test$beta_sgm22)*test$beta_sgm21[2]*test$beta_sgm22[3]

ptheta4_12 <- rep(2*test$beta_m1[6],nrow(x))

ptheta5_12 <- rep(2*test$beta_m1[6],nrow(x))

ptheta_12 <- cbind(x[,1],ptheta1_12,ptheta2_12,ptheta3_12,ptheta4_12,ptheta5_12)

thetas <- cbind(x,1/(1+pi_intm),exp(X_sig2%*%test$beta_sgm21),exp(X_sig2%*%test$beta_sgm22),X_m%*%test$beta_m1,X_m%*%test$beta_m2)

##############################################################################################################################
##############################################################################################################################
######## Simulation Study for parametric EM 
##############################################################################################################################
##############################################################################################################################
set.seed(123)
n <- 400 # set the sample size for each simulated data
nnsim <- 500 # set the number of simulated datasets

######## for p = 1  #### 
x <- runif(n)
pi_1 <- function(x){
  return(1/(1+exp(x))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(1 -0.8*x)
}

m_2 <- function(x){
  return(5 +2*x)
}

sigma_1 <- function(x){
  return(exp(-0.5+0.5*x))
}

sigma_2 <- function(x){
  return(exp(-0.7-0.2*x))
}

rnpnormmix <- function(x,p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}

p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))

X_pi <- X_m <- X_sig2 <- cbind(rep(1,n),x)
#X_m <- cbind(rep(1,n),x,x^2,x[,1]*x[,2])
#initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,2)
#initl_beta_m1 <- initl_beta_m2 <- rep(0.5,2)
#test <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
betapis <- betam1s <- betam2s <- betasig21s <- betasig22s <- NULL

for (iisim in 1:nnsim){
  y <- rnpnormmix(x,p,mu,sgm)
  initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- initl_beta_m1 <- initl_beta_m2 <- rep(0.5,2)
  out <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
  if(mean(1/(1+exp(-X_pi%*%out$beta_pi))) < 0.5) {
    betapis <- cbind(betapis,out$beta_pi)
    betam1s <- cbind(betam1s,out$beta_m1)
    betam2s <- cbind(betam2s,out$beta_m2)
    betasig21s <- cbind(betasig21s,out$beta_sgm21)
    betasig22s <- cbind(betasig22s,out$beta_sgm22)
  }
  else{
    betapis <- cbind(betapis,-out$beta_pi)
    betam1s <- cbind(betam1s,out$beta_m2)
    betam2s <- cbind(betam2s,out$beta_m1)
    betasig21s <- cbind(betasig21s,out$beta_sgm22)
    betasig22s <- cbind(betasig22s,out$beta_sgm21)
  }
}

apply(betapis,1,mean)
apply(betam1s,1,mean)
apply(betam2s,1,mean)
apply(betasig21s,1,mean)
apply(betasig22s,1,mean)

apply(betapis,1,sd)
apply(betam1s,1,sd)
apply(betam2s,1,sd)
apply(betasig21s,1,sd)
apply(betasig22s,1,sd)

##############################################################################################################################

################################### test for second order degree polynomial mean functions for p==1
set.seed(123)
n <- 400 # set the sample size for each simulated data
nnsim <- 500 # set the number of simulated datasets

######## for p = 1  #### 
x <- runif(n)
pi_1 <- function(x){
  return(1/(1+exp(x))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(1 -0.8*x-x^2)
}

m_2 <- function(x){
  return(4 +2*x+3*x^2)
}

sigma_1 <- function(x){
  return(exp(-0.5+0.5*x))
}

sigma_2 <- function(x){
  return(exp(-0.7-0.2*x))
}

rnpnormmix <- function(x,p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}

p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))

X_pi <-  X_sig2 <- cbind(rep(1,n),x)
X_m <- cbind(rep(1,n),x,x^2)
#X_m <- cbind(rep(1,n),x,x^2,x[,1]*x[,2])
#initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,2)
#initl_beta_m1 <- initl_beta_m2 <- rep(0.5,2)
#test <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
betapis <- betam1s <- betam2s <- betasig21s <- betasig22s <- NULL

for (iisim in 1:nnsim){
  y <- rnpnormmix(x,p,mu,sgm)
  initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,2) 
  initl_beta_m1 <- initl_beta_m2 <- rep(0.5,3)
  out <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
  if(mean(1/(1+exp(-X_pi%*%out$beta_pi))) < 0.5) {
    betapis <- cbind(betapis,out$beta_pi)
    betam1s <- cbind(betam1s,out$beta_m1)
    betam2s <- cbind(betam2s,out$beta_m2)
    betasig21s <- cbind(betasig21s,out$beta_sgm21)
    betasig22s <- cbind(betasig22s,out$beta_sgm22)
  }
  else{
    betapis <- cbind(betapis,-out$beta_pi)
    betam1s <- cbind(betam1s,out$beta_m2)
    betam2s <- cbind(betam2s,out$beta_m1)
    betasig21s <- cbind(betasig21s,out$beta_sgm22)
    betasig22s <- cbind(betasig22s,out$beta_sgm21)
  }
}

apply(betapis,1,mean)
apply(betam1s,1,mean)
apply(betam2s,1,mean)
apply(betasig21s,1,mean)
apply(betasig22s,1,mean)

apply(betapis,1,sd)
apply(betam1s,1,sd)
apply(betam2s,1,sd)
apply(betasig21s,1,sd)
apply(betasig22s,1,sd)

##############################################################################################################################

################################### test for linear functions for p=2

set.seed(123)
n <- 400 # set the sample size for each simulated data
nnsim <- 500 # set the number of simulated datasets

######## for p = 1  #### 
x <- cbind(runif(n),runif(n))
pi_1 <- function(x){
  return(1/(1+exp(x[,1]+2*x[,2]))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(1 -2*x[,1]+x[,2])
}

m_2 <- function(x){
  return(4 +x[,1]+3*x[,2])
}

sigma_1 <- function(x){
  return(exp(-1+0.4*x[,1]+0.8*x[,2]))
}

sigma_2 <- function(x){
  return(exp(-0.5-0.5*x[,1]-0.5*x[,2]))
}

rnpnormmix <- function(x,p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}

p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))

X_pi <- X_m <- X_sig2 <- cbind(rep(1,n),x)
#X_m <- cbind(rep(1,n),x,x^2,x[,1]*x[,2])
#initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,2)
#initl_beta_m1 <- initl_beta_m2 <- rep(0.5,2)
#test <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
betapis <- betam1s <- betam2s <- betasig21s <- betasig22s <- NULL

for (iisim in 1:nnsim){
  y <- rnpnormmix(x,p,mu,sgm)
  initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- initl_beta_m1 <- initl_beta_m2 <- rep(0.5,3)
  out <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
  if(mean(1/(1+exp(-X_pi%*%out$beta_pi))) < 0.5) {
    betapis <- cbind(betapis,out$beta_pi)
    betam1s <- cbind(betam1s,out$beta_m1)
    betam2s <- cbind(betam2s,out$beta_m2)
    betasig21s <- cbind(betasig21s,out$beta_sgm21)
    betasig22s <- cbind(betasig22s,out$beta_sgm22)
  }
  else{
    betapis <- cbind(betapis,-out$beta_pi)
    betam1s <- cbind(betam1s,out$beta_m2)
    betam2s <- cbind(betam2s,out$beta_m1)
    betasig21s <- cbind(betasig21s,out$beta_sgm22)
    betasig22s <- cbind(betasig22s,out$beta_sgm21)
  }
}

apply(betapis,1,mean)
apply(betam1s,1,mean)
apply(betam2s,1,mean)
apply(betasig21s,1,mean)
apply(betasig22s,1,mean)

apply(betapis,1,sd)
apply(betam1s,1,sd)
apply(betam2s,1,sd)
apply(betasig21s,1,sd)
apply(betasig22s,1,sd)

##############################################################################################################################

################################### test for second order degree polynomial mean functions for p==2


set.seed(123)
n <- 200 # set the sample size for each simulated data
nnsim <- 500 # set the number of simulated datasets

######## for p = 1  #### 
x <- cbind(runif(n),runif(n))
pi_1 <- function(x){
  return(1/(1+exp(x[,1]+2*x[,2]))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(1 -2*x[,1]+x[,2])
}

m_2 <- function(x){
  return(4 +x[,1]+3*x[,2])
}

sigma_1 <- function(x){
  return(exp(-0.1+0*x[,1]+0*x[,2]))
}

sigma_2 <- function(x){
  return(exp(-0.5-0*x[,1]-0*x[,2]))
}

rnpnormmix <- function(x,p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}

p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))

X_pi <-  X_sig2 <- cbind(rep(1,n),x)
X_m <- cbind(rep(1,n),x)
#X_m <- cbind(rep(1,n),x,x^2,x[,1]*x[,2])
#initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,2)
#initl_beta_m1 <- initl_beta_m2 <- rep(0.5,2)
#test <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
betapis <- betam1s <- betam2s <- betasig21s <- betasig22s <- NULL

for (iisim in 1:nnsim){
  y <- rnpnormmix(x,p,mu,sgm)
  initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,3)
  initl_beta_m1 <- initl_beta_m2 <- rep(0.5,3)
  out <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
  if(mean(1/(1+exp(-X_pi%*%out$beta_pi))) < 0.5) {
    betapis <- cbind(betapis,out$beta_pi)
    betam1s <- cbind(betam1s,out$beta_m1)
    betam2s <- cbind(betam2s,out$beta_m2)
    betasig21s <- cbind(betasig21s,out$beta_sgm21)
    betasig22s <- cbind(betasig22s,out$beta_sgm22)
  }
  else{
    betapis <- cbind(betapis,-out$beta_pi)
    betam1s <- cbind(betam1s,out$beta_m2)
    betam2s <- cbind(betam2s,out$beta_m1)
    betasig21s <- cbind(betasig21s,out$beta_sgm22)
    betasig22s <- cbind(betasig22s,out$beta_sgm21)
  }
}

apply(betapis,1,mean)
apply(betam1s,1,mean)
apply(betam2s,1,mean)
apply(betasig21s,1,mean)
apply(betasig22s,1,mean)

apply(betapis,1,sd)
apply(betam1s,1,sd)
apply(betam2s,1,sd)
apply(betasig21s,1,sd)
apply(betasig22s,1,sd)

#out <- normmixEM_contvar(X_pi,X_m,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_sgm21,initl_sgm22)
initl_sgm <- apply(sgm,2,mean)^2
initl_sgm21 <- initl_sgm[1]
initl_sgm22 <- initl_sgm[2]
betapis <- betam1s <- betam2s <- sig21s <- sig22s <- NULL

set.seed(456)
for (iisim in 1:nnsim){
  y <- rnpnormmix(x,p,mu,sgm)
  initl_beta_pi <- rep(0.5,3)
  initl_beta_m1 <- initl_beta_m2 <- rep(0.5,3)
  #out1 <- normmixEM_contvar(X_pi,X_m,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_sgm21,initl_sgm22)
  out1 <- normmixEM_contvar(X_pi,X_m,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,0.5,0.5)
  if(mean(1/(1+exp(-X_pi%*%out1$beta_pi))) < 0.5) {
    betapis <- cbind(betapis,out1$beta_pi)
    betam1s <- cbind(betam1s,out1$beta_m1)
    betam2s <- cbind(betam2s,out1$beta_m2)
    sig21s <- cbind(sig21s,out1$sgm21)
    sig22s <- cbind(sig22s,out1$sgm22)
  }
  else{
    betapis <- cbind(betapis,-out1$beta_pi)
    betam1s <- cbind(betam1s,out1$beta_m2)
    betam2s <- cbind(betam2s,out1$beta_m1)
    sig21s <- cbind(sig21s,out1$sgm22)
    sig22s <- cbind(sig22s,out1$sgm21)
  }
}

apply(betapis,1,mean)
apply(betam1s,1,mean)
apply(betam2s,1,mean)
apply(sig21s,1,mean)
apply(sig22s,1,mean)

apply(betapis,1,sd)
apply(betam1s,1,sd)
apply(betam2s,1,sd)
apply(sig21s,1,sd)
apply(sig22s,1,sd)



c11 <- c(3.8024,1.1908 ,3.3145 ,-2.1283 ,1.3053 ,-2.1231 ,0.8184 ,1.5679 ,-0.8524 ,-1.0289 ,-1.1618)

c22 <- c(4,1,3,-2,1.5,-2,0.8,1.6,-1,-1,-1)

##############################################################################################################################

################################### Compare performance using cv criteria
set.seed(123)
n <- 800 # set the sample size for each simulated data
nnsim <- 500 # set the number of simulated datasets

pi_1 <- function(x){
  return(1/(1+exp(x[,1]+2*x[,2]))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(1 -2*x[,1]+x[,2])
}

m_2 <- function(x){
  return(4 +x[,1]+3*x[,2])
}

sigma_1 <- function(x){
  return(exp(-0.1+0*x[,1]+0*x[,2]))
}

sigma_2 <- function(x){
  return(exp(-0.5-0*x[,1]-0*x[,2]))
}

rnpnormmix <- function(x,p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}

cv_error1 <- cv_error2 <- rep(0,nnsim)

for (iisim in 1:nnsim){
  
  x <- cbind(runif(n),runif(n))
  p <- cbind(pi_1(x),pi_2(x))
  mu <- cbind(m_1(x),m_2(x))
  sgm <- cbind(sigma_1(x),sigma_2(x))
  y <- rnpnormmix(x,p,mu,sgm)
  
  #### split simulated data into training and test sets 
  sample <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE, prob=c(0.7,0.3))
  x_train  <- x[sample, ]
  x_test   <- x[!sample, ]
  y_train  <- y[sample]
  y_test   <- y[!sample]
  
  X_train_pi <-  X_train_sig2 <- cbind(rep(1,nrow(x_train)),x_train)
  X_train_m <- cbind(rep(1,nrow(x_train)),x_train)
  X_test_pi <-  X_test_sig2 <- cbind(rep(1,nrow(x_test)),x_test)
  X_test_m <- cbind(rep(1,nrow(x_test)),x_test)
  
  ## set initial values 
  initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,3)
  initl_beta_m1 <- initl_beta_m2 <- rep(0.5,3)
  out <- normmixEM(X_train_pi,X_train_m,X_train_sig2,y_train,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
  out_contvar <- normmixEM_contvar(X_train_pi,X_train_m,y_train,initl_beta_pi,initl_beta_m1,initl_beta_m2,0.5,0.5)
  
  est_pi1 <- c(1/(1+exp(-X_test_pi%*%out$beta_pi)))
  est_pi2 <- 1 - est_pi1
  est_m1 <- c(X_test_m%*%out$beta_m1)
  est_m2 <- c(X_test_m%*%out$beta_m2)
  est_sd1 <- c(exp(X_test_sig2%*%out$beta_sgm21/2))
  est_sd2 <- c(exp(X_test_sig2%*%out$beta_sgm22/2))
  
  est_ctv_pi1 <- c(1/(1+exp(-X_test_pi%*%out_contvar$beta_pi)))
  est_ctv_pi2 <- 1 - est_ctv_pi1
  est_ctv_m1 <- c(X_test_m%*%out_contvar$beta_m1)
  est_ctv_m2 <- c(X_test_m%*%out_contvar$beta_m2)
  est_ctv_sd1 <- c(sqrt(out_contvar$sgm21))
  est_ctv_sd2 <- c(sqrt(out_contvar$sgm22))
  
  r1 <- est_pi1*dnorm(y_test,mean = est_m1,sd = est_sd1)/(est_pi1*dnorm(y_test,mean = est_m1,sd = est_sd1)+est_pi2*dnorm(y_test,mean = est_m2,sd = est_sd2))
  r2 <- est_pi2*dnorm(y_test,mean = est_m2,sd = est_sd2)/(est_pi1*dnorm(y_test,mean = est_m1,sd = est_sd1)+est_pi2*dnorm(y_test,mean = est_m2,sd = est_sd2))
  r1_ctv <- est_ctv_pi1*dnorm(y_test,mean = est_ctv_m1,sd = est_ctv_sd1)/(est_ctv_pi1*dnorm(y_test,mean = est_ctv_m1,sd = est_ctv_sd1)+est_ctv_pi2*dnorm(y_test,mean = est_ctv_m2,sd = est_ctv_sd2))
  r2_ctv <- est_ctv_pi2*dnorm(y_test,mean = est_ctv_m2,sd = est_ctv_sd2)/(est_ctv_pi1*dnorm(y_test,mean = est_ctv_m1,sd = est_ctv_sd1)+est_ctv_pi2*dnorm(y_test,mean = est_ctv_m2,sd = est_ctv_sd2))
  
  cv_error1[iisim] <- (y_test - r1*est_m1 - r2*est_m2)^2
  cv_error2[iisim] <- (y_test - r1_ctv*est_ctv_m1 - r2_ctv*est_ctv_m2)^2
  
}

print(mean(cv_error1))
print(mean(cv_error2))
print(sd(cv_error1))
print(sd(cv_error2))


##############################################################################################################################

################################### residual plots for normal mixture with constant variance

set.seed(123)
n <- 5000 # set the sample size for each simulated data

pi_1 <- function(x){
  return(1/(1+exp(x[,1]+2*x[,2]))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(1 -2*x[,1]+x[,2])
}

m_2 <- function(x){
  return(4 +x[,1]+3*x[,2])
}

sigma_1 <- function(x){
  return(exp(-0.1+0*x[,1]+0*x[,2]))
}

sigma_2 <- function(x){
  return(exp(-0.5-0*x[,1]-0*x[,2]))
}

rnpnormmix <- function(x,p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}

x <- cbind(runif(n),runif(n))
p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))
y <- rnpnormmix(x,p,mu,sgm)

X_pi <-  X_m <- X_sig2 <- cbind(rep(1,nrow(x)),x)

initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,3)
initl_beta_m1 <- initl_beta_m2 <- rep(0.5,3)
#out <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
out_contvar <- normmixEM_contvar(X_pi,X_m,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,0.5,0.5)
est_ctv_pi1 <- c(1/(1+exp(-X_pi%*%out_contvar$beta_pi)))
est_ctv_pi2 <- 1 -est_ctv_pi1 
est_ctv_m1 <- c(X_m%*%out_contvar$beta_m1)
est_ctv_m2 <- c(X_m%*%out_contvar$beta_m2)
est_ctv_sd1 <- c(sqrt(out_contvar$sgm21))
est_ctv_sd2 <- c(sqrt(out_contvar$sgm22))

r1_ctv <- est_ctv_pi1*dnorm(y,mean = est_ctv_m1,sd = est_ctv_sd1)/(est_ctv_pi1*dnorm(y,mean = est_ctv_m1,sd = est_ctv_sd1)+est_ctv_pi2*dnorm(y,mean = est_ctv_m2,sd = est_ctv_sd2))
r2_ctv <- est_ctv_pi2*dnorm(y,mean = est_ctv_m2,sd = est_ctv_sd2)/(est_ctv_pi1*dnorm(y,mean = est_ctv_m1,sd = est_ctv_sd1)+est_ctv_pi2*dnorm(y,mean = est_ctv_m2,sd = est_ctv_sd2))


fit_ctv <- r1_ctv*est_ctv_m1 + r2_ctv*est_ctv_m2
resid2_ctv <- (y - fit_ctv)^2
plot(fit_ctv,resid2_ctv,xlab="Predicted Value",ylab="Squared Residual",main="Normal Mixtures with Constant Variance")

##############################################################################################################################

################################### residual plots for normal mixture with functional variance (but assume constant variance) 

set.seed(123)
n <- 5000 # set the sample size for each simulated data

pi_1 <- function(x){
  return(1/(1+exp(x[,1]+2*x[,2]))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(1 -2*x[,1]+x[,2])
}

m_2 <- function(x){
  return(4 +x[,1]+3*x[,2])
}

sigma_1 <- function(x){
  return(exp(-0.5+0.4*x[,1]+0.8*x[,2]))
}

sigma_2 <- function(x){
  return(exp(-0.5*x[,1]-0.5*x[,2]))
}

rnpnormmix <- function(x,p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}

x <- cbind(runif(n),runif(n))
p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))
y <- rnpnormmix(x,p,mu,sgm)

X_pi <-  X_m <- X_sig2 <- cbind(rep(1,nrow(x)),x)

initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,3)
initl_beta_m1 <- initl_beta_m2 <- rep(0.5,3)
#out <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
out_contvar <- normmixEM_contvar(X_pi,X_m,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,0.5,0.5)
est_ctv_pi1 <- c(1/(1+exp(-X_pi%*%out_contvar$beta_pi)))
est_ctv_pi2 <- 1 -est_ctv_pi1 
est_ctv_m1 <- c(X_m%*%out_contvar$beta_m1)
est_ctv_m2 <- c(X_m%*%out_contvar$beta_m2)
est_ctv_sd1 <- c(sqrt(out_contvar$sgm21))
est_ctv_sd2 <- c(sqrt(out_contvar$sgm22))

r1_ctv <- est_ctv_pi1*dnorm(y,mean = est_ctv_m1,sd = est_ctv_sd1)/(est_ctv_pi1*dnorm(y,mean = est_ctv_m1,sd = est_ctv_sd1)+est_ctv_pi2*dnorm(y,mean = est_ctv_m2,sd = est_ctv_sd2))
r2_ctv <- est_ctv_pi2*dnorm(y,mean = est_ctv_m2,sd = est_ctv_sd2)/(est_ctv_pi1*dnorm(y,mean = est_ctv_m1,sd = est_ctv_sd1)+est_ctv_pi2*dnorm(y,mean = est_ctv_m2,sd = est_ctv_sd2))


fit_ctv <- r1_ctv*est_ctv_m1 + r2_ctv*est_ctv_m2
resid2_ctv <- (y - fit_ctv)^2
plot(fit_ctv,resid2_ctv,xlab="Predicted Value",ylab="Squared Residual",main="Normal Mixtures with Functional Variance")



################################################################################################################################################
############################################## test the variance 
#############################################################################################################################################
set.seed(345)
n <- 400 # set the sample size for each simulated data
nnsim <- 500 # set the number of simulated datasets

######## for p = 1  #### 
x <- runif(n)
pi_1 <- function(x){
  return(1/(1+exp(4*x-0.8))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(1 -0.8*x-x^2)
}

m_2 <- function(x){
  return(4 +2*x+3*x^2)
}

sigma_1 <- function(x){
  return(exp(0*x))
}

sigma_2 <- function(x){
  return(exp(-0.2-0*x))
}

mix_var <- function(x){
  return(pi_1(x)*sigma_1(x)^2 + pi_2(x)*sigma_2(x)^2 + pi_1(x)*m_1(x)^2 + pi_2(x)*m_2(x)^2 - (pi_1(x)*m_1(x)+pi_2(x)*m_2(x))^2)
}


rnpnormmix <- function(x,p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}

p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))

X_pi <-  X_sig2 <- cbind(rep(1,n),x)
X_m <- cbind(rep(1,n),x,x^2)

y <- rnpnormmix(x,p,mu,sgm)
initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,2) 
initl_beta_m1 <- initl_beta_m2 <- rep(0.5,3)
out <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
est_pi1 <- c(1/(1+exp(-X_pi%*%out$beta_pi)))
est_pi2 <- 1 -est_pi1 
est_m1 <- c(X_m%*%out$beta_m1)
est_m2 <- c(X_m%*%out$beta_m2)
est_var1 <- c(exp(X_sig2%*%out$beta_sgm21))
est_var2 <- c(exp(X_sig2%*%out$beta_sgm22))
#est_mixvar <- est_pi1*est_var1 + est_pi2*est_var2 + est_pi1*est_m1^2 + est_pi2*est_m2^2 - (est_pi1*est_m1+est_pi2*est_m2)^2
est_mixvar <- est_pi1*est_var1 + est_pi2*est_var2 


plot(x,est_mixvar,ylab="Variance",main="Variance of Y",ylim=c(0.5,2))
xp=seq(0,1,0.01)
lines(xp,mix_var2(xp),type="l",col='red')


out_contvar <- normmixEM_contvar(X_pi,X_m,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,0.5,0.5)
est_ctv_pi1 <- c(1/(1+exp(-X_pi%*%out_contvar$beta_pi)))
est_ctv_pi2 <- 1 -est_ctv_pi1 
est_ctv_m1 <- c(X_m%*%out_contvar$beta_m1)
est_ctv_m2 <- c(X_m%*%out_contvar$beta_m2)
est_ctv_var1 <- c(out_contvar$sgm21)
est_ctv_var2 <- c(out_contvar$sgm22)
#est_ctv_mixvar <- est_ctv_pi1*est_ctv_var1 + est_ctv_pi2*est_ctv_var2 + est_ctv_pi1*est_ctv_m1^2 + est_ctv_pi2*est_ctv_m2^2 - (est_ctv_pi1*est_ctv_m1+est_ctv_pi2*est_ctv_m2)^2
est_ctv_mixvar <- est_ctv_pi1*est_ctv_var1 + est_ctv_pi2*est_ctv_var2  

plot(x,est_ctv_mixvar,ylab="Variance",main="Variance of Y (constant)",ylim=c(0.5,2))
lines(xp,mix_var2(xp),type="l",col='red')
#lines(xp,mix_var2(xp)+0.2,type="l",col='red',lty='dotted')



mix_var2 <- function(x){
  return(pi_1(x)*sigma_1(x)^2 + pi_2(x)*sigma_2(x)^2)
}
plot(xp,mix_var(xp),type='l')

pi_1 <- function(x){
  return(1/(1+exp(4*x-0.8))) 
}
mix_var3 <- function(x){
  return(pi_1(x)*1 + pi_2(x)*0.67032)
}
plot(xp,mix_var3(xp),type='l')
plot(xp,pi_1(xp),type='l')

xpp <- seq(0,10,0.01)
plot(xpp,mix_var2(xpp),type='l')
plot(x,y)
lines(xp,m_1(xp),type='l')
lines(xp,m_2(xp),type='l',col='red')

##############################################################################################################################

################################### test for second order degree polynomial mean functions for p=1


set.seed(234)
n <- 400 # set the sample size for each simulated data
nnsim <- 500 # set the number of simulated datasets

######## for p = 1  #### 
x <- cbind(runif(n))
pi_1 <- function(x){
  return(1/(1+exp(4*x-0.8))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(1 -0.8*x-x^2)
}

m_2 <- function(x){
  return(5 +2*x+3*x^2)
}

sigma_1 <- function(x){
  return(exp(0.2+0*x))
}

sigma_2 <- function(x){
  return(exp(-0.2-0*x))
}

mix_var <- function(x){
  return(pi_1(x)*sigma_1(x)^2 + pi_2(x)*sigma_2(x)^2 + pi_1(x)*m_1(x)^2 + pi_2(x)*m_2(x)^2 - (pi_1(x)*m_1(x)+pi_2(x)*m_2(x))^2)
}

mix_var2 <- function(x){
  return(pi_1(x)*sigma_1(x)^2 + pi_2(x)*sigma_2(x)^2)
}

true_betapi <- c(0.8,-4)
true_betam1 <- c(1,-0.8,-1)
true_betam2 <- c(5,2,3)
true_betasig21 <- c(0.4,0)
true_betasig22 <- c(-0.4,0)

rnpnormmix <- function(x,p,m,sgm){
  y <- NULL
  for (i in 1:nrow(p)){
    y[i] <- rnormmix(1,p[i,],m[i,],sgm[i,])
  }
  return(y)
}

p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm <- cbind(sigma_1(x),sigma_2(x))

X_pi <-  X_sig2 <- cbind(rep(1,n),x)
X_m <- cbind(rep(1,n),x,x^2)
#X_m <- cbind(rep(1,n),x,x^2,x[,1]*x[,2])
#initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,2)
#initl_beta_m1 <- initl_beta_m2 <- rep(0.5,2)
#test <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
betapis <- betam1s <- betam2s <- betasig21s <- betasig22s <- NULL

for (iisim in 1:nnsim){
  y <- rnpnormmix(x,p,mu,sgm)
  initl_beta_pi <-  initl_beta_sgm21 <- initl_beta_sgm22 <- rep(0.5,2)
  initl_beta_m1 <- initl_beta_m2 <- rep(0.5,3)
  out <- normmixEM(X_pi,X_m,X_sig2,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
  if(out$beta_pi[2] < 0) {
    betapis <- cbind(betapis,out$beta_pi)
    betam1s <- cbind(betam1s,out$beta_m1)
    betam2s <- cbind(betam2s,out$beta_m2)
    betasig21s <- cbind(betasig21s,out$beta_sgm21)
    betasig22s <- cbind(betasig22s,out$beta_sgm22)
  }
  else{
    betapis <- cbind(betapis,-out$beta_pi)
    betam1s <- cbind(betam1s,out$beta_m2)
    betam2s <- cbind(betam2s,out$beta_m1)
    betasig21s <- cbind(betasig21s,out$beta_sgm22)
    betasig22s <- cbind(betasig22s,out$beta_sgm21)
  }
}

apply(betapis,1,mean)
apply(betam1s,1,mean)
apply(betam2s,1,mean)
apply(betasig21s,1,mean)
apply(betasig22s,1,mean)

apply(betapis,1,sd)
apply(betam1s,1,sd)
apply(betam2s,1,sd)
apply(betasig21s,1,sd)
apply(betasig22s,1,sd)

xp <- seq(0,1,0.01)
Xp <- cbind(rep(1,length(xp)),xp)
Xpm <- cbind(rep(1,length(xp)),xp,xp^2)
estpis <-1/(1+exp(t(-Xp%*%betapis)))

pi1_hat <- function(x){
  return(1/(1+exp(4.083*x-0.827))) 
}


plot(xp,pi_1(xp),type='l',col='red',ylim=c(0,1),xlab='x',ylab='y')
#lines(xp,pi1_hat(xp))
lines(xp,apply(estpis,2,mean),type='l')
lines(xp,apply(estpis, 2, function(x) quantile(x, probs = .975)),type='l',lty='dotted')
lines(xp,apply(estpis, 2, function(x) quantile(x, probs = .025)),type='l',lty='dotted')
legend(x = "topright",          # Position
       legend = c("estimated proportion function", "true proportion function","95% confidence intervals"),  # Legend texts
       lty = c(1, 1, 3),           # Line types
       col = c(1, 2, 1),           # Line colors
       lwd = 2)
title(main="Estimated proportion function, true proportion function and confidence interval",cex.main = 0.5)

estm1s <- t(Xpm%*%betam1s)
estm2s <- t(Xpm%*%betam2s)
plot(xp,m_1(xp),type='l',col='red',ylim=c(-3,15),xlab='x',ylab='y')
lines(xp,m_2(xp),type='l',col='red')
lines(xp,apply(estm1s,2,mean),type='l')
lines(xp,apply(estm1s, 2, function(x) quantile(x, probs = .975)),type='l',lty='dotted')
lines(xp,apply(estm1s, 2, function(x) quantile(x, probs = .025)),type='l',lty='dotted')
lines(xp,apply(estm2s,2,mean),type='l')
lines(xp,apply(estm2s, 2, function(x) quantile(x, probs = .975)),type='l',lty='dotted')
lines(xp,apply(estm2s, 2, function(x) quantile(x, probs = .025)),type='l',lty='dotted')
legend(x = "topleft",          # Position
       legend = c("estimated mean functions", "true mean functions","95% confidence intervals"),  # Legend texts
       lty = c(1, 1, 3),           # Line types
       col = c(1, 2, 1),           # Line colors
       lwd = 2)
title(main="Estimated mean functions, true mean functions and confidence interval",cex.main = 0.5)

estsig21s <- exp(t(Xp%*%betasig21s))
estsig22s <- exp(t(Xp%*%betasig22s))
plot(xp,sigma_1(xp)^2,type='l',col='red',ylim=c(0,12),xlab='x',ylab='y')
lines(xp,sigma_2(xp)^2,type='l',col='red')
lines(xp,apply(estsig21s,2,mean),type='l')
lines(xp,apply(estsig21s, 2, function(x) quantile(x, probs = .975)),type='l',lty='dotted')
lines(xp,apply(estsig21s, 2, function(x) quantile(x, probs = .025)),type='l',lty='dotted')
lines(xp,apply(estsig22s,2,mean),type='l')
lines(xp,apply(estsig22s, 2, function(x) quantile(x, probs = .975)),type='l',lty='dotted')
lines(xp,apply(estsig22s, 2, function(x) quantile(x, probs = .025)),type='l',lty='dotted')
legend(x = "topleft",          # Position
       legend = c("estimated variance functions", "true variance functions","95% confidence intervals"),  # Legend texts
       lty = c(1, 1, 3),           # Line types
       col = c(1, 2, 1),           # Line colors
       lwd = 2)
title(main="Estimated variance functions, true variance functions and confidence interval",cex.main = 0.5)

estmixvar <- estpis*estsig21s + (1-estpis)*estsig22s
plot(xp,mix_var2(xp),type='l',col='red',xlab='x',ylab='Weighted variance',ylim=c(0.5,1.8))
lines(xp,apply(estmixvar,2,mean),type='l')
lines(xp,apply(estmixvar, 2, function(x) quantile(x, probs = .975)),type='l',lty='dotted')
lines(xp,apply(estmixvar, 2, function(x) quantile(x, probs = .025)),type='l',lty='dotted')
legend(x = "topright",          # Position
       legend = c("estimated weighted variance", "true weighted variance","95% confidence intervals"),  # Legend texts
       lty = c(1, 1, 3),           # Line types
       col = c(1, 2, 1),           # Line colors
       lwd = 2)
title(main="Estimated weighted, true weighted variance and confidence interval
      (M2, setting 2) ",cex.main = 0.5)


mse_betapi <- rowSums((betapis - true_betapi)^2)/nnsim 
mse_betam1 <- rowSums((betam1s - true_betam1)^2)/nnsim 
mse_betam2 <- rowSums((betam2s - true_betam2)^2)/nnsim 
mse_betasig21 <- rowSums((betasig21s - true_betasig21)^2)/nnsim
mse_betasig22 <- rowSums((betasig22s - true_betasig22)^2)/nnsim

#out <- normmixEM_contvar(X_pi,X_m,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_sgm21,initl_sgm22)
initl_sgm <- apply(sgm,2,mean)^2
initl_sgm21 <- initl_sgm[1]
initl_sgm22 <- initl_sgm[2] 
betapis <- betam1s <- betam2s <- sig21s <- sig22s <- NULL

set.seed(234)
x <- cbind(runif(n))
for (iisim in 1:nnsim){
  y <- rnpnormmix(x,p,mu,sgm)
  initl_beta_pi <- rep(0.5,2)
  initl_beta_m1 <- initl_beta_m2 <- rep(0.5,3)
  #out1 <- normmixEM_contvar(X_pi,X_m,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_sgm21,initl_sgm22)
  out1 <- normmixEM_contvar(X_pi,X_m,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,0.5,0.5)
  if(out1$beta_pi[2] < 0) {
    betapis <- cbind(betapis,out1$beta_pi)
    betam1s <- cbind(betam1s,out1$beta_m1)
    betam2s <- cbind(betam2s,out1$beta_m2)
    sig21s <- cbind(sig21s,out1$sgm21)
    sig22s <- cbind(sig22s,out1$sgm22)
  }
  else{
    betapis <- cbind(betapis,-out1$beta_pi)
    betam1s <- cbind(betam1s,out1$beta_m2)
    betam2s <- cbind(betam2s,out1$beta_m1)
    sig21s <- cbind(sig21s,out1$sgm22)
    sig22s <- cbind(sig22s,out1$sgm21)
  }
}
betasig21s <- rbind(log(sig21s),rep(0,length(sig21s)))
betasig22s <- rbind(log(sig22s),rep(0,length(sig22s)))

apply(betapis,1,mean)
apply(betam1s,1,mean)
apply(betam2s,1,mean)
apply(sig21s,1,mean)
apply(sig22s,1,mean)

apply(betapis,1,sd)
apply(betam1s,1,sd)
apply(betam2s,1,sd)
apply(sig21s,1,sd)
apply(sig22s,1,sd)


################################################################################################################################################################
############################################   Real data application                               ###########################################
################################################################################################################################################################
ydata <- as.numeric(CO2GDP2005$`2005 [YR2005]`[ CO2GDP2005$`Series Name`=="CO2 emissions (metric tons per capita)"])/5
#xdata <- as.numeric(CO2GDP2005$`2005 [YR2005]`[ CO2GDP2005$`Series Name`=="GDP per capita (current US$)"])/10000
xdata <- log(as.numeric(CO2GDP2005$`2005 [YR2005]`[ CO2GDP2005$`Series Name`=="GDP per capita (current US$)"]))

plot(xdata,ydata,ylab='CO2 per capita (in metric tons)',xlab='log GDP per capita',main="Scatter plot of the CO2-GDP data")

#### use the traditional mixture regression to find the initial values 
library(mixreg)
test_mixreg <- mixreg(xdata,ydata,ncomp =2)
test_mixreg$parmat

ndata <- length(ydata)
Xdata <- cbind(rep(1,ndata),xdata)
# set the initial data 
initl_beta_pi <-  c(0.3,-0.5) 
initl_beta_m1 <- c(-4.323719,0.6881191)
initl_beta_m2 <- c(-1.556115,0.2560951)
initl_beta_sgm21 <- c(-0.5,0.1)
initl_beta_sgm22 <- c(-3.5,0.1)

output <- normmixEM(Xdata,Xdata,Xdata,ydata,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
print(output)

##### find the bootstrap 95% CI for betas
est_pi1 <- c(1/(1+exp(-Xdata%*%output$beta_pi)))
est_pi2 <- 1 -est_pi1 
est_m1 <- c(Xdata%*%output$beta_m1)
est_m2 <- c(Xdata%*%output$beta_m2)
est_var1 <- c(exp(Xdata%*%output$beta_sgm21))
est_var2 <- c(exp(Xdata%*%output$beta_sgm22))
p <- cbind(est_pi1,est_pi2)
mu <- cbind(est_m1 ,est_m2)
sgm <- cbind(sqrt(est_var1),sqrt(est_var2))

betapis <- betam1s <- betam2s <- betasig21s <- betasig22s <- NULL

set.seed(3668)
for (iisim in 1:ndata){
  y <- rnpnormmix(Xdata,p,mu,sgm)
  initl_beta_pi <-  output$beta_pi  
  initl_beta_m1 <- output$beta_m1
  initl_beta_m2 <- output$beta_m2
  initl_beta_sgm21 <- output$beta_sgm21
  initl_beta_sgm22 <- output$beta_sgm22
  out <- normmixEM(Xdata,Xdata,Xdata,y,initl_beta_pi,initl_beta_m1,initl_beta_m2,initl_beta_sgm21,initl_beta_sgm22)
    betapis <- cbind(betapis,out$beta_pi)
    betam1s <- cbind(betam1s,out$beta_m1)
    betam2s <- cbind(betam2s,out$beta_m2)
    betasig21s <- cbind(betasig21s,out$beta_sgm21)
    betasig22s <- cbind(betasig22s,out$beta_sgm22)
}

print(apply(betapis,1,function(x) quantile(x, probs = .975)))
print(apply(betapis,1,function(x) quantile(x, probs = .025)))
print(apply(betam1s,1,function(x) quantile(x, probs = .975)))
print(apply(betam1s,1,function(x) quantile(x, probs = .025)))
print(apply(betam2s,1,function(x) quantile(x, probs = .975)))
print(apply(betam2s,1,function(x) quantile(x, probs = .025)))
print(apply(betasig21s,1,function(x) quantile(x, probs = .975)))
print(apply(betasig21s,1,function(x) quantile(x, probs = .025)))
print(apply(betasig22s,1,function(x) quantile(x, probs = .975)))
print(apply(betasig22s,1,function(x) quantile(x, probs = .025)))













