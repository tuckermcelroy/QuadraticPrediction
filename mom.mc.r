mom.mc <- function(sims,horz)
{

	T <- dim(sims)[1]
	T2 <- floor(T/2)
	mc <- dim(sims)[2]

	lags <- seq(0,horz-1)
	lags <- c(-rev(lags),lags[-1])
	
	# 1st moment
	mu <- mean(sims[1,])

	# 2nd moment
	mom2 <- array(0,length(lags))
	for(h in lags) { mom2[h+horz] <- mean((sims[T2,]-mu)*(sims[T2+h,]-mu)) }

	# 3rd moment
	mom3 <- array(0,rep(length(lags),2))
	for(h1 in lags)
	{
	for(h2 in lags)
	{
		mom3[h1+horz,h2+horz] <- 
			mean((sims[T2,]-mu)*(sims[T2+h1,]-mu)*(sims[T2+h2,]-mu))
	}}

	# 4th moment
	mom4 <- array(0,rep(length(lags),3))
	for(h1 in lags)
	{
	for(h2 in lags)
	{
	for(h3 in lags)
	{
		mom4[h1+horz,h2+horz,h3+horz] <- 
			mean((sims[T2,]-mu)*(sims[T2+h1,]-mu)*(sims[T2+h2,]-mu)*(sims[T2+h3,]-mu))
	}}}

	return(list(mu,mom2,mom3,mom4))
}

