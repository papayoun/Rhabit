
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
#' @param lag_inter integer specifying which points of cov_array will be used
#' for the interpolation, default to 10
#' xgrid and ygrid, of dimensions (length(xgrid),length(ygrid),length(beta)).
#' @export
simLangevinMM <- function(beta, gamma2 = 1, times, loc0,
                          cov_list = NULL, grad_fun = NULL) {
  checkCovGrad(cov_list, grad_fun)
  nb_obs <- length(times)
  xy <- matrix(loc0, nb_obs, 2)
  dt <- diff(times)
  for (t in 2:nb_obs) {
    cat("\rSimulating Langevin process...", round(100 * t / nb_obs), "%")
    grad <- gradLogUD(beta = beta, loc = xy[t - 1, ], cov_list = cov_list,
                      grad_fun = grad_fun, check = F)
    rand_part <- stats::rnorm(2, 0, sqrt(gamma2 * dt[t - 1]))
    xy[t, ] <- xy[t - 1, ] + 0.5 * grad * gamma2 * dt[t - 1] + rand_part
  }
  cat("\n")
  return(data.frame(x = xy[, 1], y = xy[, 2], t = times))
}
