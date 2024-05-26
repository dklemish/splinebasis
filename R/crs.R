crs <- function(x, knots, incl_intercept=FALSE){
  #' Set up a cubic regression spline basis for a single variable at user 
  #' specified knot locations.
  #' 
  #' @param x A vector of data.
  #' @param knots A vector of knot locations defining the cubic regression spline basis.
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
  h <- diff(knots)

  # Define B matrix per Table 5.1 of Wood 2ed
  B <- matrix(0, nrow=K-2, ncol=K-2)
  for(i in 1:(K-2)){
    B[i,i] <- (h[i] + h[i+1])/3
  }
  
  for(i in 1:(K-3)){
    B[i, i+1] <- h[i+1]/6
    B[i+1, i] <- h[i+1]/6
  }
  
  # Define D matrix per Table 5.1 of Wood 2ed
  D <- matrix(0, nrow=K-2, ncol=K)
  for(i in 1:(K-2)){
    D[i,i]    <- 1/h[i]
    D[i, i+1] <- -1/h[i] - 1/h[i+1]
    D[i, i+2] <- 1/h[i+1]
  }
  
  # Define F matrix
  Fminus <- solve(B) %*% D
  Fmat   <- rbind(rep(0, K), Fminus, rep(0, K))
  
  # Define a functions
  n <- length(x)
  
  a_minus <- matrix(0, nrow=n, ncol=K-1)
  a_plus  <- matrix(0, nrow=n, ncol=K-1)
  for(j in 1:(K-1)){
    a_minus[, j] <- (knots[j+1] - x)/h[j]
    a_plus[, j]  <- (x - knots[j])/h[j]
  }
  
  # Define c functions
  c_minus <- matrix(0, nrow=length(x), ncol=K-1)
  c_plus  <- matrix(0, nrow=length(x), ncol=K-1)
  for(j in 1:(K-1)){
    c_minus[, j] <- ((knots[j+1] - x)^3/h[j] - h[j]*(knots[j+1] - x)) / 6
    c_plus[, j]  <- ((x - knots[j])^3/h[j] - h[j]*(x - knots[j]) ) / 6
  }
  
  # Define implicit basis functions b
  b <- matrix(0, nrow=n, ncol=K)
  
  # j indexes spline functions, i x-location, ki for knot interval
  for(j in 1:K){
    for(ki in 1:(K-1)){
      for(i in 1:n){
        if((x[i] >= knots[ki]) & (x[i] <= knots[ki+1])){
          if(ki == j){
            b[i, j] <- a_minus[i, ki] + c_minus[i, ki]*Fmat[ki, j] + c_plus[i, ki]*Fmat[ki+1, j]          
          } else if(ki==(j-1)){
            b[i, j] <- a_plus[i, ki]  + c_minus[i, ki]*Fmat[ki, j] + c_plus[i, ki]*Fmat[ki+1, j]          
          } else{
            b[i, j] <- c_minus[i, ki]*Fmat[ki, j] + c_plus[i, ki]*Fmat[ki+1, j]
          }
        }    
      }
    }
  }
  
  if(incl_intercept){
    N <- cbind(1, x, b[,-1])
  } else{
    N <- cbind(x, b)
  }
  
  return(N)
}