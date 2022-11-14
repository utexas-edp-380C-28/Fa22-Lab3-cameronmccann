#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 3: An Overview of Linear Regression 
#        and Data Generation
#
#' This function generates a matrix of multivariate data based on the 
#' sample size (n), mean column vector (mu), and covariance matrix (Sigma)
#' 
#' @param n Sample size (or number of rows) 
#' @param mu A column vector of means 
#' @param Sigma Covariance matrix 
#' @return Returns a matrix of the generated data 
rmvnorm <- function(n, mu, Sigma) {
  # Error message if dimensions do not match
  try(if(nrow(mu) != nrow(Sigma))
    stop("Input dimensions do not match"))
  
  # Populate Z
  Z <- matrix(1, n, nrow(mu))
  
  # Assign standard normal deviates in Z
  for (i in 1:length(mu)) {
    Z[, i] <- rnorm(n)
  }
  
  # Generate data
  generated_data <- matrix(1, n, 1) %*% t(mu) + Z %*% chol(Sigma)
  
  return(generated_data)
}

