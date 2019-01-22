
#' Simulate from the Langevin movement model
#'
#' This function is based on the Euler approximation.
#'
#' @name simLangevinMM
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
#' @param keep_grad should gradient values at simulated points be kept?
#' This avoids to use covGradAtLocs later on
#' xgrid and ygrid, of dimensions (length(xgrid),length(ygrid),length(beta)).
#' @export
simLangevinMM <- function(beta, gamma2 = 1, times, loc0,
                          cov_list = NULL, grad_fun = NULL, silent = F,
                          keep_grad = F) {
  checkCovGrad(cov_list, grad_fun)
  nb_obs <- length(times)
  
  # Matrix of simulated locations
  xy <- matrix(NA, nb_obs, 2)
  xy[1,] <- loc0
  
  # Vector of time intervals
  dt <- diff(times)
  
  # Number of covariates
  J <- length(beta)
  
  # Initialise gradient array
  grad_array <- NULL
  if (keep_grad){
    grad_array <- matrix(ncol = 2 * J, nrow = nb_obs)
    colnames(grad_array) <- paste(rep(paste0("grad_c", 1:J), rep(2, J)),
                                  rep(c("x", "y"), J), sep = "_")
  }
  
  # Simulate locations
  for (t in 2:nb_obs) {
    if (!silent)
      cat("\rSimulating Langevin process...", round(100 * t / nb_obs), "%")
    
    # Only keep subpart of covariate maps (to save memory)
    cov_list_tmp <- lapply(cov_list, getGridZoom, x0 = xy[t - 1, ])
    
    # Compute covariate gradients at current location
    grad_val <- bilinearGrad(loc = xy[t - 1, ], cov_list = cov_list_tmp)
    
    if (keep_grad)
      grad_array[t - 1, ] <- as.numeric(grad_val)
    
    # Gradient of utilisation distribution
    grad <- grad_val %*% beta
    
    # Simulate next location
    rand_part <- stats::rnorm(2, 0, sqrt(gamma2 * dt[t - 1]))
    xy[t, ] <- xy[t - 1, ] + 0.5 * grad * gamma2 * dt[t - 1] + rand_part
  }
  
  if (!silent)
    cat("\n")
  
  # Data frame of simulated locations
  main_df <- data.frame(x = xy[, 1], y = xy[, 2], t = times)
  
  if (keep_grad){
    cov_list_tmp <- lapply(cov_list, getGridZoom, x0 = xy[nb_obs - 1, ])
    
    grad_val <- bilinearGrad(loc = xy[t - 1, ], cov_list = cov_list_tmp)
    
    grad_array[nb_obs, ] <- as.numeric(grad_val)
    main_df <- cbind.data.frame(main_df, grad_array)
  }
  
  return(main_df)
}
