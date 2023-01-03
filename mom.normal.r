mom.normal <- function(gamma)
{
	horz <- (length(gamma)-1)/2
	lags <- seq(0,horz-1)
	lags <- c(-rev(lags),lags[-1])
	
	# 1st moment
	mu <- 0

	# 2nd moment
	mom2 <- array(0,length(lags))
	for(h in lags) { mom2[h+horz] <- gamma[abs(h)+1] }

	# 3rd moment
	mom3 <- array(0,rep(length(lags),2))
 
	# 4th moment
	mom4 <- array(0,rep(length(lags),3))
 
	return(list(mu,mom2,mom3,mom4))
}




