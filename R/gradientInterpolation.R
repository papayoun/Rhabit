#' Interpolate 2D covariate
#'
#' At the moment, this is based on the function \code{interp.surface} from
#' the package fields (bilinear interpolation).
#'
#' @param locs Point where the covariate should be interpolated
#' @param xgrid Grid on which the covariate is known
#' @param ygrid Grid on which the covariate is known
#' @param cov_mat Matrix of values of the covariate at the points given by
#' xgrid and ygrid.
#' @return Interpolated value of the covariate at the point xy.
#' @details Note that covmat needs to rotated (as e.g. with "image"), so you
#' might need to use something like
#'   covmat <- t(apply(as.matrix(covraster),2,rev))
#' before passing it to this function
interpCov <- function(locs, xgrid, ygrid, cov_mat) {
  locs <- matrix(locs, ncol = 2)
  cov_xyz <- list(x = xgrid, y = ygrid, z = cov_mat)
  return(fields::interp.surface(cov_xyz, locs))
}

#' Gradient of the log of the utilisation distribution
#'
#' @param beta Vector of resource selection coefficients
#' @param loc Point at which the gradient should be evaluated (2d vector)
#' @param grad_fun Optional list of functions taking a 2d vector and returning
#'  a two 2d vector for the gradient
#' @param xgrid Grid on which the covariates are known
#' @param ygrid Grid on which the covariates are known
#' @param cov_array Array of values of the covariates at the points given by
#' xgrid and ygrid, of dimensions (length(xgrid),length(ygrid),length(beta)).
#'
#' @return Gradient of the log-UD in loc.
gradLogUD <- function(beta, loc, xgrid, ygrid, cov_array, grad_fun = NULL) {
  J <- length(beta)
  grad_vals <- sapply(1:J, function(j)
    if(is.null(grad_fun[[j]])){
      cov_mat <- cov_array[,,j]
      return(nloptr::nl.grad(fn = interpCov, x0 = loc, xgrid = xgrid,
                      ygrid = ygrid, cov_mat = cov_mat))
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
#' @param xgrid Grid on which the covariates are known
#' @param ygrid Grid on which the covariates are known
#' @param cov_array Array of values of the covariates at the points given by
#' xgrid and ygrid, of dimensions (length(xgrid),length(ygrid),length(beta)).
#' @param lag_inter integer specifying which points of cov_array will be used
#' for the interpolation, default to 10
#' @param grad_fun Optional list of functions taking a 2d vector and returning
#' a two 2d vector for the gradient
#' @return Three-dimensional array of gradients of covariate fields.
#' The rows index time, the columns are the dimensions (x and y),
#' and the layers index the covariates.
#' @export
covGradAtLocs <- function(locs, xgrid, ygrid, cov_array,
                          lag_inter = 10, grad_fun = NULL) {
  J <- dim(cov_array)[3]
  x_size <- length(xgrid)
  y_size <- length(ygrid)
  grad_array <- do.call(function(...) abind::abind(..., along = 3),
                        lapply(1:dim(cov_array)[3], function(j){
                          if(is.null(grad_fun[[j]])){
                            return(t(apply(locs, 1, function(x){
                              x_pos <- findInterval(x[1], xgrid)
                              y_pos <- findInterval(x[2], ygrid)
                              min_indx <- max(1, x_pos - lag_inter)
                              max_indx <- min(x_size, x_pos + lag_inter)
                              x_inter <- min_indx:max_indx
                              min_indy <- max(1, y_pos - lag_inter)
                              max_indy <- min(y_size, y_pos + lag_inter)
                              y_inter <- min_indy:max_indy
                              nloptr::nl.grad(fn = interpCov, x0 = x,
                                              xgrid = xgrid[x_inter],
                                              ygrid = ygrid[y_inter],
                                              cov_mat = cov_array[x_inter,
                                                                  y_inter, j])
                            })))
                          }
                          else
                            return(apply(locs, 1, grad_fun[[j]]))
                        }))
  return(grad_array)
}
