
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
#' @param grad_fun Optional list of functions taking a 2d vector and returning
#'  a two 2d vector for the gradient
#' @param xgrid Grid on which the covariates are known
#' @param ygrid Grid on which the covariates are known
#' @param cov_array Array of values of the covariates at the points given by
#' @param lag_inter integer specifying which points of cov_array will be used
#' for the interpolation, default to 10
#' xgrid and ygrid, of dimensions (length(xgrid),length(ygrid),length(beta)).
#' @export
simLangevinMM <- function(beta, gamma2 = 1, times, loc0,
                          xgrid, ygrid, cov_array,
                          lag_inter = 10, grad_fun = NULL) {
  nb_obs <- length(times)
  xy <- matrix(loc0, nb_obs, 2)
  x_size <- length(xgrid)
  y_size <- length(ygrid)
  dt <- diff(times)
  for (t in 2:nb_obs) {
    cat("\rSimulating Langevin process...", round(100 * t / nb_obs), "%")
    x_pos <- findInterval(xy[t - 1, 1], xgrid)
    y_pos <- findInterval(xy[t - 1, 2], ygrid)
    x_inter <- max(1, x_pos - lag_inter):min(x_size, x_pos + lag_inter)
    y_inter <- max(1, y_pos - lag_inter):min(y_size, y_pos + lag_inter)
    grad <- gradLogUD(beta = beta, loc = xy[t - 1, ], xgrid = xgrid[x_inter],
                      ygrid = ygrid[y_inter],
                      cov_array = cov_array[x_inter, y_inter, ],
                      grad_fun = grad_fun)
    rand_part <- stats::rnorm(2, 0, sqrt(gamma2 * dt[t - 1]))
    xy[t, ] <- xy[t - 1, ] + 0.5 * grad * gamma2 * dt[t - 1] + rand_part
  }
  cat("\n")
  return(data.frame(x = xy[, 1], y = xy[, 2], t = times))
}
