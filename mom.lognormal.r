mom.lognormal <- function(gamma)
{

	horz <- (length(gamma)-1)/2
	lags <- seq(0,horz-1)
	lags <- c(-rev(lags),lags[-1])
	
	# 1st moment
	mu <- exp(gamma[1]/2)

	# 2nd moment
	mom2 <- array(0,length(lags))
	for(h in lags) { mom2[h+horz] <- mu^2*(exp(gamma[abs(h)+1])-1) }

	# 3rd moment
	mom3 <- array(0,rep(length(lags),2))
	for(h1 in lags)
	{
	for(h2 in lags)
	{
		mom3[h1+horz,h2+horz] <- mu^3*(exp(gamma[abs(h1)+1]+gamma[abs(h2)+1]+gamma[abs(h1-h2)+1]) -
			exp(gamma[abs(h1)+1]) - exp(gamma[abs(h2)+1]) - exp(gamma[abs(h1-h2)+1]) + 2)
	}}

	# 4th moment
	mom4 <- array(0,rep(length(lags),3))
	for(h1 in lags)
	{
	for(h2 in lags)
	{
	for(h3 in lags)
	{
		mom4[h1+horz,h2+horz,h3+horz] <- mu^4*(exp(gamma[abs(h1)+1]+gamma[abs(h2)+1]+gamma[abs(h3)+1]+
			gamma[abs(h1-h2)+1]+gamma[abs(h1-h3)+1]+gamma[abs(h2-h3)+1]) - 
			exp(gamma[abs(h1-h2)+1]+gamma[abs(h1-h3)+1]+gamma[abs(h2-h3)+1]) -
			exp(gamma[abs(h2)+1]+gamma[abs(h3)+1]+gamma[abs(h2-h3)+1]) -
			exp(gamma[abs(h1)+1]+gamma[abs(h3)+1]+gamma[abs(h1-h3)+1]) -
			exp(gamma[abs(h1)+1]+gamma[abs(h2)+1]+gamma[abs(h1-h2)+1]) +
			exp(gamma[abs(h2-h3)+1]) + exp(gamma[abs(h1-h3)+1]) + exp(gamma[abs(h1-h2)+1]) +
			exp(gamma[abs(h3)+1]) + exp(gamma[abs(h2)+1]) + exp(gamma[abs(h1)+1]) -3)
	}}}

	return(list(mu,mom2,mom3,mom4))
}



