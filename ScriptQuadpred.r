## Script for Quadratic Prediction Project

setwd("C:\\Users\\neide\\OneDrive\\Documents\\Research\\NonlinearPred")

#source("ARMAauto.r")
#source("getGCD.r")
source("quad.pred.r")
source("quad.cast.r")
source("mom.normal.r")
source("mom.lognormal.r")
source("mom.sqrtnormal.r")
source("mom.nonparam.r")
source("mom.mc.r")
#source("garch.sim.r")


##################################################################
##  nonparametric estimation of automoments 

T <- 100

# example 1: Gaussian WN
x <- rnorm(T)

# example 2: quadratic WN process
z <- rnorm(T+1)
x <- z[1:T]*z[2:(T+1)]

horz <- 5
lags <- seq(0,horz-1)
lags <- c(-rev(lags),lags[-1])

moms <- mom.nonparam(x,horz)

#####################################################
## computation of automoments for an asymmetric GARCH(1,1) process

monte <- 10^3
alpha <- c(.3,.4)
x.sims <- NULL
for(i in 1:monte)
{
  x.sim <- garch.sim(alpha,n=100,g=rexp.ctr)
  x.sims <- cbind(x.sims,x.sim)
}

out <- mom.mc(x.sims,5)

###########################################################
## computation of automoments for a lognormal process

T <- 5
H <- 1
horz <- T+H
lags <- seq(0,horz-1)
lags <- c(-rev(lags),lags[-1])

# process 1: an AR(1)
phi <- -.8
sigma <- 1
gamma <- ARMAauto(phi,NULL,(2*horz))*sigma^2

# process 2: an MA(1)
theta <- .9
sigma <- 1
gamma <- ARMAauto(NULL,theta,(2*horz))*sigma^2

# process 3: an AR(2)
rho <- .8
phi1 <- sqrt(3)*rho
phi2 <- -rho^2
sigma <- .5
gamma <- ARMAauto(c(phi1,phi2),NULL,(2*horz))*sigma^2

moms <- mom.lognormal(gamma)


###########################################################
## computation of automoments for a gaussian-chisquare process
# higher values of alpha weight Gaussian part more, alpha=0 is pure chi-square

T <- 5
H <- 1
horz <- T+H
lags <- seq(0,horz-1)
lags <- c(-rev(lags),lags[-1])

# process 1: an AR(1)
alpha <- 0
phi <- -.9
sigma <- 1
gamma <- ARMAauto(phi,NULL,(2*horz))*sigma^2
gamma <- gamma/gamma[1]

# process 2: an MA(1)
theta <- -1
alpha <- sqrt(abs(theta))
sigma <- 1
gamma <- ARMAauto(NULL,theta,(2*horz))*sigma^2
gamma <- gamma/gamma[1]

# process 3: an AR(2)
alpha <- 2
rho <- .99
omega <- pi/6
phi1 <- 2*cos(omega)*rho
phi2 <- -rho^2
sigma <- .1
gamma <- ARMAauto(c(phi1,phi2),NULL,(2*horz))*sigma^2
gamma <- gamma/gamma[1]

# process 4: seasonal MA
alpha <- 0
theta <- rep(1,3)
sigma <- 1
gamma <- ARMAauto(NULL,theta,(2*horz))*sigma^2
gamma <- gamma/gamma[1]

moms <- mom.sqrtnormal(gamma,alpha)


######################################################
##  quadratic prediction for h-step ahead forecasting
#	requires mom2, mom3, mom4 up to lag T+H-1

output <- quad.pred(T,H,moms)
mses <- output[[1]]
quad.coef <- output[[2]]
lin.coef <- output[[3]]
quad.mat <- output[[4]]
print(mses[2]/mses[1])


## construct forecast 

# lognormal case simulation
z.lognormal <- t(chol(quad.mat)) %*% rnorm(T)
x.lognormal <- exp(z.lognormal) - mu
data <- x.lognormal

# quadratic hermite case simulation
z.hermite <- t(chol(quad.mat)) %*% rnorm(T)
x.hermite <- alpha*z.hermite + (z.hermite^2 - 1)
data <- x.hermite

# forecast and plots
forecast <- quad.cast(data,T,H,moms)
plot(ts(c(data,forecast),start=1990,frequency=1),col=2,lty=1,ylab="",xlab="",lwd=2)
lines(ts(c(data,NA),start=1990,frequency=1),lwd=2)


################################################
### Data Analysis: wolfer sunspots

wolfer <- read.table("wolfer.dat")
wolfer.raw <- ts(wolfer,start=1749,frequency=12)
wolfer <- sqrt(wolfer.raw)

plot(wolfer)

n <- length(wolfer)
mu <- mean(wolfer)
sig <- sd(wolfer)
gamma.hat <- acf(wolfer,lag=n-1,type="covariance")$acf[,,1]
kappa.hat <- pacf(wolfer,lag=n-1)$acf[,,1]

### parametric
# AR id
alpha <- .05
crit <- qnorm(1-alpha/2)
K.n <- 1 + floor(3*sqrt(log(n,base=10)))
kappa.test <- kappa.hat/sqrt(log(n,base=10)/n)
k <- 1
while(k < (n-K.n))
{
  if(max(abs(kappa.test[k:(k+K.n-1)])) < crit) 
  { 
    p.hat <- k-1
    k <- n-K.n
  } else { k <- k+1 }
}
print(p.hat)

phi.ar <- solve(toeplitz(gamma.hat[1:p.hat])) %*% gamma.hat[2:(p.hat+1)]
sig2 <- gamma.hat[1] - t(gamma.hat[2:(p.hat+1)]) %*%
  solve(toeplitz(gamma.hat[1:p.hat])) %*% gamma.hat[2:(p.hat+1)]

ar.poly <- c(1,-1*phi.ar)
ar.roots <- polyroot(ar.poly)
Mod(ar.roots)
2*pi/Arg(ar.roots)

T <- p.hat
H <- 1
horz <- T+H
gamma <- ARMAauto(phi.ar,NULL,(2*horz))*sig2[1,1]
rho <- gamma/gamma[1]

moms <- mom.sqrtnormal(rho,2*mu/sig)

### nonparametric

T <- 30
H <- 1
horz <- T+H
moms <- mom.nonparam(wolfer.raw,horz)


## get results

output <- quad.pred(T,H,moms)
mses <- output[[1]]
quad.coef <- output[[2]]
lin.coef <- output[[3]]
quad.mat <- output[[4]]
print(mses[2]/mses[1])

# forecast and plots
forecast <- quad.cast(wolfer.raw,T,H,moms)
plot(ts(c(wolfer.raw,forecast),start=1749,frequency=12),col=2,lty=1,ylab="",xlab="",lwd=2)
lines(ts(c(wolfer.raw,NA),start=1749,frequency=12),lwd=2)



