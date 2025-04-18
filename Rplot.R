################################################################################################################
##### Plots of rejection proportions 
###############################################################################################################
par(mfrow = c(2, 2))

beta <- c(0,0.25,0.5,0.75,1)
p1 <- c(0.06,0.092,0.276,0.47,0.708)
p2 <- c(0.088,0.116,0.246,0.41,0.662)

plot(beta, p1,type = "b", frame = TRUE,pch = 18,
     col = "red", xlab = expression(beta), ylab = "Proportion of Rejection",ylim = c(0,1),main = expression(paste("Comparison of Rejection Proportions,",alpha,"=0.05,","n=100")))
# Add a second line
lines(beta, p2, col = "blue",pch = 18, lty = 2,type = "b")
# Add a legend to the plot
legend(0.01,0.99, legend=c("Model 1", "Model 2"),
       col=c("red", "blue"), lty = 1:2,cex=0.8,bty = "n")

p1 <- c(0.104,0.164,0.38,0.602,0.816)
p2 <- c(0.138,0.192,0.33,0.59,0.766)

plot(beta, p1,type = "b", frame = TRUE,pch = 18,
     col = "red", xlab = expression(beta), ylab = "Proportion of Rejection",ylim = c(0,1),main = expression(paste("Comparison of Rejection Proportions,",alpha,"=0.10,","n=100")))
# Add a second line
lines(beta, p2, col = "blue",pch = 18, lty = 2,type = "b")
# Add a legend to the plot
legend(0.01,0.99, legend=c("Model 1", "Model 2"),
       col=c("red", "blue"), lty = 1:2,cex=0.8,bty = "n")

p1 <- c(0.06,0.164,0.48,0.85,0.974)
p2 <- c(0.066,0.156,0.402,0.742,0.934)

plot(beta, p1,type = "b", frame = TRUE,pch = 18,
     col = "red", xlab = expression(beta), ylab = "Proportion of Rejection",ylim = c(0,1),main = expression(paste("Comparison of Rejection Proportions,",alpha,"=0.05,","n=200")))
# Add a second line
lines(beta, p2, col = "blue",pch = 18, lty = 2,type = "b")
# Add a legend to the plot
legend(0.01,0.99, legend=c("Model 1", "Model 2"),
       col=c("red", "blue"), lty = 1:2,cex=0.8,bty = "n")

p1 <- c(0.11,0.24,0.586,0.896,0.986)
p2 <- c(0.112,0.234,0.52,0.81,0.958)

plot(beta, p1,type = "b", frame = TRUE,pch = 18,
     col = "red", xlab = expression(beta), ylab = "Proportion of Rejection",ylim = c(0,1),main = expression(paste("Comparison of Rejection Proportions,",alpha,"=0.10,","n=200")))
# Add a second line
lines(beta, p2, col = "blue",pch = 18, lty = 2,type = "b")
# Add a legend to the plot
legend(0.01,0.99, legend=c("Model 1", "Model 2"),
       col=c("red", "blue"), lty = 1:2,cex=0.8,bty = "n")
text(-0.1,0.1,"a",cex=2)

par(mfrow = c(1, 1))

################################################################################################################
##### Scatter plots of simulated samples under the Null
###############################################################################################################
library(mixtools)

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

n <- 200
set.seed(599)
x <- runif(n)
p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm1 <- cbind(rep(0.7,n),rep(0.4,n))
sgm2 <- cbind(sigma_1(x),rep(0.4,n))
y1 <- rnpnormmix(p,mu,sgm1)
y2 <- rnpnormmix(p,mu,sgm2)

par(mfrow = c(1, 2))
plot(x,y1,ylab = "y",main="Example of Simulated Dataset, Model 1")
plot(x,y2,ylab = "y",main="Example of Simulated Dataset, Model 2")
par(mfrow = c(1, 1))

################################################################################################################
##### Scatter plots of simulated samples under heteroscedasticity
###############################################################################################################

pi_1 <- function(x){
  return(1/(1+exp(0.8-2*x))) 
}

pi_2 <- function(x){
  return(1-pi_1(x))
}

m_1 <- function(x){
  return(2 -2*x )   
}

m_2 <- function(x){
  return(6 +3*x )
}

sigma_1_a1 <- function(x){
  return(0.1+0.5*exp(-1.2+1.5*x))
}
sigma_1_b1 <- function(x){
  return(0.1+0.5*(x-0.5)^2)
}

sigma_1_c1 <- function(x){
  return(0.1+0.5*abs(sin(2*pi*x)))
}
sigma_1_d1 <- function(x){
  return(0.6+0.5*(x<0.7&x>0.3))
}

sigma_1_a2 <- function(x){
  return(0.1+1*exp(-1.2+1.5*x))
}
sigma_1_b2 <- function(x){
  return(0.1+1*(x-0.5)^2)
}

sigma_1_c2 <- function(x){
  return(0.1+1*abs(sin(2*pi*x)))
}
sigma_1_d2 <- function(x){
  return(0.6+1*(x<0.7&x>0.3))
}

sigma_2 <- 0.4
n <- 200 
set.seed(599)
x <- runif(n)
p <- cbind(pi_1(x),pi_2(x))
mu <- cbind(m_1(x),m_2(x))
sgm_a1 <- cbind(sigma_1_a1(x),rep(sigma_2,n))
sgm_b1 <- cbind(sigma_1_b1(x),rep(sigma_2,n))
sgm_c1 <- cbind(sigma_1_c1(x),rep(sigma_2,n))
sgm_d1 <- cbind(sigma_1_d1(x),rep(sigma_2,n))
sgm_a2 <- cbind(sigma_1_a2(x),rep(sigma_2,n))
sgm_b2 <- cbind(sigma_1_b2(x),rep(sigma_2,n))
sgm_c2 <- cbind(sigma_1_c2(x),rep(sigma_2,n))
sgm_d2 <- cbind(sigma_1_d2(x),rep(sigma_2,n))

y_a1 <- rnpnormmix(p,mu,sgm_a1)
y_b1 <- rnpnormmix(p,mu,sgm_b1)
y_c1 <- rnpnormmix(p,mu,sgm_c1)
y_d1 <- rnpnormmix(p,mu,sgm_d1)
y_a2 <- rnpnormmix(p,mu,sgm_a2)
y_b2 <- rnpnormmix(p,mu,sgm_b2)
y_c2 <- rnpnormmix(p,mu,sgm_c2)
y_d2 <- rnpnormmix(p,mu,sgm_d2)

par(mfrow = c(2, 2))
plot(x,y_a1,ylab = "y",main = expression(paste("Example of Simulated Dataset,","H"[1]^a,",",lambda,"=0.5")))
plot(x,y_b1,ylab = "y",main = expression(paste("Example of Simulated Dataset,","H"[1]^b,",",lambda,"=0.5")))
plot(x,y_c1,ylab = "y",main = expression(paste("Example of Simulated Dataset,","H"[1]^c,",",lambda,"=0.5")))
plot(x,y_d1,ylab = "y",main = expression(paste("Example of Simulated Dataset,","H"[1]^d,",",lambda,"=0.5")))
par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
plot(x,y_a2,ylab = "y",main = expression(paste("Example of Simulated Dataset,","H"[1]^a,",",lambda,"=1")))
plot(x,y_b2,ylab = "y",main = expression(paste("Example of Simulated Dataset,","H"[1]^a,",",lambda,"=1")))
plot(x,y_c2,ylab = "y",main = expression(paste("Example of Simulated Dataset,","H"[1]^a,",",lambda,"=1")))
plot(x,y_d2,ylab = "y",main = expression(paste("Example of Simulated Dataset,","H"[1]^a,",",lambda,"=1")))
par(mfrow = c(1, 1))




