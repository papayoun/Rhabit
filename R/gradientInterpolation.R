#' Interpolate 2D covariate
#'
#' At the moment, this is based on the function \code{interp.surface} from
#' the package fields (bilinear interpolation).
#'
#' @param locs Point where the covariate should be interpolated
#' @param x_grid Grid on which the covariate is known
#' @param y_grid Grid on which the covariate is known
#' @param cov_mat Matrix of values of the covariate at the points given by
#' x_grid and y_grid.
#' @return Interpolated value of the covariate at the point xy.
#' @details Note that covmat needs to rotated (as e.g. with "image"), so you
#' might need to use something like
#'   covmat <- t(apply(as.matrix(covraster),2,rev))
#' before passing it to this function
interpCov <- function(locs, x_grid, y_grid, cov_mat) {
  locs <- matrix(locs, ncol = 2)
  cov_xyz <- list(x = x_grid, y = y_grid, z = cov_mat)
  return(fields::interp.surface(cov_xyz, locs))
}

#' Gradient of the log of the utilisation distribution
#'
#' @param beta Vector of resource selection coefficients
#' @param loc Point at which the gradient should be evaluated (2d vector)
#' @param cov_list List of J (number of covariates) "raster like" elements.
#' A raster like element is a 3 elements list with named elements
#' 1) "x" a vector of increasing x locations (at which the covariate is sampled)
#' 2) "y" a vector of increasing y locations (at which the covariate is sampled)
#' 3) "z" a size(x)*size(y) matrix giving covariate values at location (x, y)
#' @param grad_fun Optional list of functions taking a 2d vector and returning
#' a two 2d vector for the gradient
#' @return Gradient of the log-UD in loc.
gradLogUD <- function(beta, loc, cov_list = NULL, grad_fun = NULL, check = T) {
  if(check){
    checkCovGrad(cov_list, grad_fun)
  }
  J <- length(beta)
  grad_vals <- sapply(1:J, function(j)
    if(is.null(grad_fun[[j]])){
      x_grid <- cov_list[[j]]$x
      y_grid <- cov_list[[j]]$y
      cov_mat <- cov_list[[j]]$z
      return(nloptr::nl.grad(fn = interpCov, x0 = loc, x_grid = x_grid,
                      y_grid = y_grid, cov_mat = cov_mat))
    }
    else
      return(grad_fun[[j]](loc))
    )
  return(grad_vals %*% beta)
}

#' Gradient of covariate field
#'
#' @name covGradAtLocs
#' @param locs Matrix of locations where the gradient should be evaluated
#' @param cov_list List of J (number of covariates) "raster like" elements.
#' A raster like element is a 3 elements list with named elements
#' 1) "x" a vector of increasing x locations (at which the covariate is sampled)
#' 2) "y" a vector of increasing y locations (at which the covariate is sampled)
#' 3) "z" a size(x)*size(y) matrix giving covariate values at location (x, y)
#' @param lag_inter integer specifying which points of cov_array will be used
#' for the interpolation, default to 10
#' @param grad_fun Optional list of functions taking a 2d vector and returning
#' a two 2d vector for the gradient
#' @return Three-dimensional array of gradients of covariate fields.
#' The rows index time, the columns are the dimensions (x and y),
#' and the layers index the covariates.
#' @export
covGradAtLocs <- function(locs, cov_list = NULL, grad_fun = NULL) {
  checkCovGrad(cov_list, grad_fun)
  J <- length(cov_list)
  grad_array <- do.call(function(...) abind::abind(..., along = 3),
                        lapply(1:J, function(j){
                          if(is.null(grad_fun[[j]])){
                            x_grid <- cov_list[[j]]$x
                            y_grid <- cov_list[[j]]$y
                            cov_mat <- cov_list[[j]]$z
                            return(t(apply(locs, 1, function(x){
                              nloptr::nl.grad(fn = interpCov, x0 = x,
                                              x_grid = x_grid,
                                              y_grid = y_grid,
                                              cov_mat = cov_mat)
                            })))
                          }
                          else
                            return(apply(locs, 1, grad_fun[[j]]))
                        }))
  return(grad_array)
}
