tcs <- function(x, knots, incl_intercept=FALSE){
  #' Set up a truncated cubic spline basis for a single variable at user 
  #' specified knot locations.
  #' 
  #' @param x A vector of data.
  #' @param knots A vector of knot locations defining the cubic spline basis.
  #' @param incl_intercept A boolean indicating whether the resulting spline 
  #' basis functions should include a column of all 1's. While a column of 1's 
  #' is normally included in the natural cubic spline basis, it can be useful
  #' to exclude it if the resulting spline basis functions will be used in 
  #' a regression formula as in lm or glm.
  #' @return A matrix of K columns where K is the number of knots (or K-1 if 
  #' an intercept is excluded)
  #' @details This implpementation is based on section 5.2.1 of 'Elements of 
  #' Statistical Learning, 2nd ed. by Hastie, Tibshirani & Friedman.
  
  K <- length(knots)
  N <- matrix(0, nrow=length(x), ncol=K+4)
  N[,1] <- 1
  N[,2] <- x
  N[,3] <- x^2
  N[,4] <- x^3
  
  for(i in 1:K){
    N[,i+4] <- pmax(0, (x-knots[i])^3) 
  }
  
  if(incl_intercept){
    N <- N[, -1]
  }
  
  return(N)
}