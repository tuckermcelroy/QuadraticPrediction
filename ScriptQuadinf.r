## Script for "Prediction for Second Order Forecastable Processes"

setwd("C:\\Users\\neide\\OneDrive\\Documents\\GitHub\\QuadraticPrediction")

source("ARMAauto.r")
source("quad.pred.r")
source("quad.cast.r")
source("mom.normal.r")
source("mom.lognormal.r")
source("mom.sqrtnormal.r")
source("mom.nonparam.r")
source("mom.mc.r")
#source("garch.sim.r")
source("specFactmvar.r")
source("polymult.r")

## Compute quadratic filter for signal+noise process

# stdev of Gaussian white noise
tau <- 1  
# moments of centered degree 1 chi square 
mu2 <- 2
mu3 <- 8
mu4 <- 48
# moving average parameter
theta <- .5

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
N <- 6
xAcf <- array(0,c(N,2,N))
Gamma <- toeplitz(c((phi^2 + (1+phi^2)^2)*sigma2^2, phi*(1+phi^2)*sigma2^2, rep(0,N-4)))
gamma0 <- c((1+phi^2)*sigma2, (1+theta^2)*mu3, rep(0,N-2))
gamma1 <- c((1+theta^2)*mu3, 2*(1+phi^2)^2*sigma2^2, 2*phi*(1+phi^2)^2*sigma2^2, rep(0,N-3))
Gamma <- rbind(gamma1[-c(1,2)],Gamma)
Gamma <- rbind(gamma0[-c(1,2)],Gamma)
Gamma <- cbind(gamma1,Gamma)
Gamma <- cbind(gamma0,Gamma)
xAcf[,1,] <- Gamma
Gamma <- toeplitz(c(phi^2*sigma2^2,phi*(1+phi^2)*sigma2^2,phi^2*sigma2^2,rep(0,N-5)))
Gamma[lower.tri(Gamma)] <- 0
gamma0 <- c(phi*sigma2,theta*mu3,theta*mu3,rep(0,N-3))
gamma1 <- c(theta*mu3,2*phi^2*sigma2^2+mu4*theta^2,phi*(1+phi^2)*sigma2^2,phi^2*sigma2^2,rep(0,N-4))
Gamma <- rbind(gamma1[-c(1,2)],Gamma)
Gamma <- rbind(gamma0[-c(1,2)],Gamma)
Gamma <- cbind(gamma1,Gamma)
Gamma <- cbind(c(gamma0[c(1,2)],rep(0,N-2)),Gamma)
xAcf[,2,] <- Gamma

