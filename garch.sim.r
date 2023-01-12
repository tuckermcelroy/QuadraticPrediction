garch.sim <- function(alpha, n=1000, b=100 ,g = rnorm, ...){
  if(n <= 0)
    stop("Sample Size must be a positive integer")

  total.n = n + b
  z = g(total.n, ...)
  z = scale(z)
  x = numeric(total.n)

  x[1:10]=z[1:10]
  for(i in 11:total.n){
    x[i] = z[i]*sqrt(alpha[1] + alpha[2]*(x[i-1]^2))
  }

  y= x[c((b+1):total.n)]
  return(y)
}

