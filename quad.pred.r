quad.pred <- function(P,H,moms)
{

  # quad.pred by Tucker McElroy (January 2023)
  #
  # Computes formulas for quadratic prediction
  #
  # Inputs:
  #   P: Number of past variables utilized for prediction
  #   H: horizon for forecasts
  #   moms: list object for autocumulants
  #     moms[[1]] <- mean
  #     moms[[2]] <- 1-array of autocovariances
  #     moms[[3]] <- 2-array of 3rd order autocumulants
  #     moms[[4]] <- 3-array of 4th order autocumulants
  # Outputs:
  #   lin.mse, quad.mse, total.mse: MSE contributions from linear and quadratic
  #   quad.coef: coefficients for quadratic predictors
  #   lin.coef: coefficients for linear predictors
  #   quad.mat: P x P covariance matrix
  
	horz <- P+H
	quad.mom <- array(0,c(P,P))
	cubic.mom <- array(0,c(P,P,P))
	quart.mom <- array(0,c(P,P,P,P))
	ylin.mom <- array(0,P)
	yquad.mom <- array(0,c(P,P))

	# indices: i,k are row and column within a block, which is indexed by j,l

	# quadratic moment
	for(i in 1:P)
	{
	for(k in 1:P)
	{
		quad.mom[i,k] <- moms[[2]][i-k+horz]
	}}

	# cubic moment
	for(i in 1:P)
	{
	for(k in 1:P)
	{
	for(l in 1:P)
	{
		cubic.mom[i,k,l] <- moms[[3]][i-l+horz,k-l+horz]
	}}}

	# quartic moment
	for(i in 1:P)
	{
	for(j in 1:P)
	{
	for(k in 1:P)
	{
	for(l in 1:P)
	{
		quart.mom[i,j,k,l] <- moms[[4]][i-l+horz,j-l+horz,k-l+horz] -
			moms[[2]][i-j+horz]*moms[[2]][k-l+horz]
	}}}}

	for(k in 1:P)
	{
		ylin.mom[k] <- moms[[2]][P+H-k+horz]
	}

	for(k in 1:P)
	{
	for(l in 1:P)
	{
		yquad.mom[k,l] <- moms[[3]][P+H-l+horz,k-l+horz]
	}}

	ind.mat <- matrix(seq(1,P^2),P,P)
	ind.omits <- ind.mat[upper.tri(ind.mat)]
	quad.mat <- quad.mom
	cubic.mat <- matrix(cubic.mom,c(P,P^2))
	cubic.mat <- cubic.mat[,-ind.omits]
	quart.mat <- matrix(quart.mom,c(P^2,P^2))
	quart.mat <- quart.mat[-ind.omits,-ind.omits]
	yquad.mat <- matrix(yquad.mom[lower.tri(ind.mat,diag=TRUE)],ncol=1)

	## check invertibility of solution
	soln.mat <- rbind(cbind(quad.mat,cubic.mat),cbind(t(cubic.mat),quart.mat))
	eigs <- eigen(soln.mat)$values
	if(min(eigs) < 10^(-5)) print(c("Singularity",min(eigs)))

	## construct solution
 	schur.mat <- quart.mat - t(cubic.mat) %*% solve(quad.mat) %*% cubic.mat
	lin.predict <- solve(quad.mat) %*% matrix(ylin.mom,ncol=1)
	quad.predict <- yquad.mat - t(cubic.mat) %*% solve(quad.mat) %*% matrix(ylin.mom,ncol=1)
	quad.coef <- solve(schur.mat) %*% quad.predict
	lin.coef <- lin.predict - solve(quad.mat) %*% cubic.mat %*% quad.coef

	## construct forecast MSE
	lin.mse <- moms[[2]][horz] - t(matrix(ylin.mom,ncol=1)) %*% lin.predict
	quad.mse <- t(quad.predict) %*% quad.coef
	total.mse <- lin.mse - quad.mse

	return(list(c(lin.mse,quad.mse,total.mse),quad.coef,lin.coef,quad.mat))

}
