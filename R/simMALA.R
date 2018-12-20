
#' Evaluate the log-RSF based on interpolated covariates
#' (used in simMALA)
#' 
#' @param locs Matrix of locations where the log-RSF should be evaluated
#' @param beta Vector of resource selection coefficients
#' @param cov_list List of J (number of covariates) "raster like" elements.
#' A raster like element is a 3 elements list with named elements
#' 1) "x" a vector of increasing x locations (at which the covariate is sampled)
#' 2) "y" a vector of increasing y locations (at which the covariate is sampled)
#' 3) "z" a size(x)*size(y) matrix giving covariate values at location (x, y)
logRSFinterp <- function(locs, beta, cov_list) {
  c <- sapply(cov_list, function(cov) {
    x_grid <- cov$x
    y_grid <- cov$y
    cov_mat <- cov$z
    interpCov(locs = locs, x_grid = x_grid, y_grid = y_grid, 
              cov_mat = cov_mat) 
  })
  return(c %*% beta)
}

#' Simulate from Metropolis-adjusted Langevin algorithm
#'
#' @name simMALA
#' 
#' @param beta Vector of resource selection coefficients
#' @param gamma2 Scalar speed parameter
#' @param times Vector of times of observations. Should have a small
#' time step between each, such that the Euler scheme is valid
#' @param loc0 Vector of coordinates (x,y) of initial location
#' @param cov_list List of J (number of covariates) "raster like" elements.
#' A raster like element is a 3 elements list with named elements
#' 1) "x" a vector of increasing x locations (at which the covariate is sampled)
#' 2) "y" a vector of increasing y locations (at which the covariate is sampled)
#' 3) "z" a size(x)*size(y) matrix giving covariate values at location (x, y)
#' @param grad_fun Optional list of functions taking a 2d vector and returning
#' a two 2d vector for the gradient
#' @param silent logical, should simulation advancement be shown?
#' 
#' @importFrom stats dnorm runif
simMALA <- function(beta, gamma2 = 1, times, loc0, cov_list = NULL, 
                    grad_fun = NULL, silent = FALSE) 
{
  checkCovGrad(cov_list, grad_fun)
  nb_obs <- length(times)
  J <- length(beta)
  
  # Initialise matrix of locations
  xy <- matrix(NA, nb_obs, 2)
  xy[1,] <- loc0
  # Time intervals
  dt <- diff(times)
  
  # Log-RSF and gradient of log-UD at initial position
  covgrad <- gradLogUD(beta = beta, loc=xy[1,], cov_list = cov_list,
                       grad_fun = grad_fun, check = FALSE)
  g <- covgrad %*% beta
  logRSF <- logRSFinterp(locs = xy[1,], beta = beta, cov_list = cov_list)
  
  # Count acceptances and rejections
  acc <- 0
  rej <- 0
  
  # Simulate from Metropolis-adjusted Langevin process
  for (t in 2:nb_obs) {
    if (!silent)
      cat("\rSimulating Langevin process...", round(100 * t / nb_obs), "%")
    
    # Sample proposed location
    rand_part <- stats::rnorm(2, 0, sqrt(gamma2 * dt[t - 1]))
    xyprime <- xy[t - 1, ] + 0.5 * g * gamma2 * dt[t - 1] + rand_part
    
    # Log-RSF and gradient of log-UD at proposed location
    covgrad <- gradLogUD(beta = beta, loc = xyprime, cov_list = cov_list, 
                         grad_fun = grad_fun, check = FALSE)
    gprime <- covgrad %*% beta
    logRSFprime <- logRSFinterp(locs = xyprime, beta = beta, cov_list = cov_list)
    
    # Log-proposals in both directions, for the acceptance ratio
    logProp <- sum(dnorm(xy[t-1,], xyprime + gamma2*dt[t-1]*gprime/2, 
                         sqrt(gamma2*dt[t-1]), log=TRUE))
    logPropPrime <- sum(dnorm(xyprime, xy[t-1,] + gamma2*dt[t-1]*g/2, 
                              sqrt(gamma2*dt[t-1]), log=TRUE))
    
    # Log acceptance ratio
    logAR <- logRSFprime + logProp - logRSF - logPropPrime
    if(log(runif(1))<logAR) {
      # If accepted, update location, gradient, and log-RSF
      xy[t,] <- xyprime
      g <- gprime
      logRSF <- logRSFprime
      acc <- acc+1
    } else {
      xy[t,] <- xy[t-1,]
      rej <- rej + 1
    } 
  }
  if (!silent)
    cat("\n")
  
  # Simulated data set
  main_df <- data.frame(x = xy[, 1], y = xy[, 2], t = times)
  
  return(list(data = main_df, acc = acc, rej = rej))
}
