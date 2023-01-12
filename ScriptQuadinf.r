## Script for "Prediction for Second Order Forecastable Processes"

library(xtable)

setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\QuadraticPrediction")

source("specFactmvar.r")
source("polymult.r")

## Compute quadratic filter for signal+noise process

# stdev of Gaussian white noise
#tau <- 0
#tau <- 1
tau <- 2
# moments of centered degree 1 chi square 
mu2 <- 2
mu3 <- 8
mu4 <- 48
# moving average parameter
theta <- -.8

## compute new parameters
rho1 <- mu2*theta/(mu2*(1+theta^2)+tau^2)
phi.pos <- (1 + sqrt(1-4*rho1^2))/(2*rho1)
phi.neg <- (1 - sqrt(1-4*rho1^2))/(2*rho1)
if(abs(phi.pos)< 1) { phi <- phi.pos } else { phi <- phi.neg }
sigma2 <- mu2*theta/phi

## check 
mu2*polymult(c(1,theta),rev(c(1,theta))) + c(0,tau^2,0)
sigma2*polymult(c(1,phi),rev(c(1,phi)))

## Compute acvf
# set truncation N > 4
N <- 10
xAcf <- array(0,c(N,2,N))
Gamma <- toeplitz(c((1+phi^2)^2*sigma2^2, 
                    phi*(1+phi^2)*sigma2^2, rep(0,N-5)))
Gamma <- cbind(c(phi*(1+phi^2)*sigma2^2,rep(0,N-4)),Gamma)
Gamma <- cbind(rep(0,N-3),Gamma)
Gamma <- cbind(rep(0,N-3),Gamma)
gamma0 <- c((1+phi^2)*sigma2, (1+theta^3)*mu3, theta^2*mu3, rep(0,N-3))
gamma1 <- c((1+theta^3)*mu3, 2*(1+phi^2)^2*sigma2^2 + mu4*(1+theta^4), 
            2*phi*(1+phi^2)*sigma2^2 + mu4*theta^3, rep(0,N-3))
gamma2 <- c(theta^2*mu3,2*phi*(1+phi^2)*sigma2^2 + mu4*theta^3,
            (phi^2 + (1+phi^2)^2)*sigma2^2 + mu4*theta^2,
            phi*(1+phi^2)*sigma2^2,rep(0,N-4))
Gamma <- rbind(gamma2,Gamma)
Gamma <- rbind(gamma1,Gamma)
Gamma <- rbind(gamma0,Gamma)
xAcf[,1,] <- Gamma
Gamma <- toeplitz(c(phi^2*sigma2^2,phi*(1+phi^2)*sigma2^2,phi^2*sigma2^2,rep(0,N-6)))
Gamma[lower.tri(Gamma)] <- 0
Gamma <- cbind(rep(0,N-3),Gamma)
Gamma <- cbind(c(phi^2*sigma2^2,rep(0,N-4)),Gamma)
Gamma <- cbind(rep(0,N-3),Gamma)
gamma0 <- c(phi*sigma2,theta*mu3,rep(0,N-2))
gamma1 <- c(theta^2*mu3,2*phi^2*sigma2^2+mu4*theta^2,phi*(1+phi^2)*sigma2^2,
            phi^2*sigma2^2,rep(0,N-4))
gamma2 <- c(theta*mu3,phi*(1+phi^2)*sigma2^2+mu4*theta,phi^2*sigma2^2,
            phi*(1+phi^2)*sigma2^2,phi^2*sigma2^2,rep(0,N-5))
Gamma <- rbind(gamma2,Gamma)
Gamma <- rbind(gamma1,Gamma)
Gamma <- rbind(gamma0,Gamma)
xAcf[,2,] <- Gamma

## compute factorization
factor <- specFactmvar(xAcf)
theta.vma <- factor[[1]][,,2:1]
sigma.vma <- factor[[2]]

## quadratic mse divided by linear mse
print(sigma.vma[1,1]/sigma2)

## filter coefficients
trunc <- 20
coeff <- theta.vma[,,2]
quad.filter <- matrix(0,nrow=trunc+1,ncol=N)
quad.filter[1,] <- coeff[1,]
for(j in 1:trunc)
{
  coeff <- -1*theta.vma[,,2] %*% coeff
  quad.filter[j+1,] <- coeff[1,]
}
rownames(quad.filter) <- seq(0,trunc)
print(xtable(quad.filter,digits=4))
      

