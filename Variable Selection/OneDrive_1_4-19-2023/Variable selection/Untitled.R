soft_threshold <- function(rho, lamda) {
  # Soft threshold function used for normalized data and lasso regression
  if (rho < -lamda) {
    return(rho + lamda)
  } else if (rho > lamda) {
    return(rho - lamda)
  } else {
    return(0)
  }
}

coordinate_descent_lasso <- function(theta, X, y, lamda = 0.01, num_iters = 100, intercept = FALSE) {
  # Coordinate gradient descent for lasso regression - for normalized data.
  # The intercept parameter allows specifying whether or not we regularize theta_0
  
  # Initialization of useful values
  m <- nrow(X)
  n <- ncol(X)
  X <- X / matrix(rep(sqrt(colSums(X^2))), nrow = m, ncol = n, byrow = FALSE)
  
  # Looping until max number of iterations
  for (i in 1:num_iters) {
    
    # Looping through each coordinate
    for (j in 1:n) {
      
      # Vectorized implementation
      X_j <- X[, j]
      y_pred <- X %*% theta
      rho <- t(X_j) %*% (y - y_pred + theta[j] * X_j)
      
      # Checking intercept parameter
      if (intercept) {
        if (j == 1) {
          theta[j] <- rho
        } else {
          theta[j] <- soft_threshold(rho, lamda)
        }
      }
      
      if (!intercept) {
        theta[j] <- soft_threshold(rho, lamda)
      }
    }
  }
  return(theta)
}

beta <- c(rnorm(10,5,5), rep(0,90))
# 0.76080677  0.80029371 -0.02458476 -0.62178759 -0.44165875 -0.50871856
# -0.68176433  0.53259264 -0.98431832  1.13879438 
X <- matrix(rnorm(10000), ncol = 100)
y = X %*% beta + rnorm(100)

beta_init <- rep(0,100)

result <- coordinate_descent_lasso(theta = beta_init,X = X, y = y,num_iters = 1000,intercept = FALSE )

round(result, 2)












