mom.nonparam <- function(data,horz)
{

lag.nlts <- function(data,lag)
{
	T <- length(data)
	if(lag == 0) { y <- data }
	if(lag > 0) { y <- c(rep(0,lag),data[1:(T-lag)]) }
	if(lag < 0) { y <- c(data[(-lag+1):T],rep(0,-lag)) }
	return(y)
}

mom.nlts <- function(data,lags)
{
	k <- length(lags)
	T <- length(data)
	data.lagged <- data
	for(i in 1:k)
	{
		data.lagged <- data.lagged * lag.nlts(data,lags[i])
	}
	mom.est <- sum(data.lagged)/T
	return(mom.est)
}

	lags <- seq(0,horz-1)
	lags <- c(-rev(lags),lags[-1])
	
	# 1st moment
	mu <- mean(data)

	# 2nd moment
	mom2 <- array(0,length(lags))
	for(h in lags) { mom2[h+horz] <- mom.nlts(data-mu,h) }

	# 3rd moment
	mom3 <- array(0,rep(length(lags),2))
	for(h1 in lags)
	{
	for(h2 in lags)
	{
		mom3[h1+horz,h2+horz] <- mom.nlts(data-mu,c(h1,h2))
	}}

	# 4th moment
	mom4 <- array(0,rep(length(lags),3))
	for(h1 in lags)
	{
	for(h2 in lags)
	{
	for(h3 in lags)
	{
		mom4[h1+horz,h2+horz,h3+horz] <- mom.nlts(data-mu,c(h1,h2,h3))
	}}}

	return(list(mu,mom2,mom3,mom4))
}



