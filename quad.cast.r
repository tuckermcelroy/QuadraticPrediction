quad.cast <- function(data,P,H,moms)
{

  # quad.cast by Tucker McElroy (January 2023)
  #
  # Applies quadratic predictors to data
  #
  # Inputs:
  #   data: time series sample of length n [requires n >= P]
  #   P: Number of past variables utilized for prediction
  #   H: horizon for forecasts
  #   moms: list object for autocumulants
  #     moms[[1]] <- mean
  #     moms[[2]] <- 1-array of autocovariances
  #     moms[[3]] <- 2-array of 3rd order autocumulants
  #     moms[[4]] <- 3-array of 4th order autocumulants
  # Outputs:
  #   cast: application of quadratic prediction formulas to
  #     most recent P observations.
  # Dependencies: quad.pred
  
	n <- length(data)
	horz <- P+H
	ind.mat <- matrix(seq(1,P^2),P,P)
	ind.omits <- ind.mat[upper.tri(ind.mat)]

	output <- quad.pred(P,H,moms)
	quad.coef <- output[[2]]
	lin.coef <- output[[3]]
	quad.mat <- output[[4]]
	
	data.now <- data[(n-P+1):n]
	cast <- t(lin.coef) %*% data.now + 
		t(quad.coef) %*% matrix((data.now %x% t(data.now) - quad.mat)[lower.tri(ind.mat,diag=TRUE)],ncol=1)

	return(cast)
}


