#' Simulate random covariate field
#'
#' @name simSpatialCov
#' @param nu smoothness parameter
#' @param rho range parameter
#' @param sigma2 variance parameter
#' @param lim Vector of limits of the map (xmin, xmax, ymin, ymax)
#' @param resol Grid resolution (defaults to 1)
#' @param mean_function function giving mean of a spatial coordinate. Must
#' be coded to take a 2 (or 3)d vector and return a scalar. By default set
#' to return 0 for every coordinates
#' @param raster_like logical specifying the return output. If TRUE, then
#' a list with xgrid, ygrid and  matrix is returned convenient for image.plot.
#' If FALSE a data.frame is returned, which might be more convenient for ggplot.
#' Default to FALSE.
#'
#' @return either a data.frame (convenient for ggplot2) or a raster like list
#' @export
simSpatialCov <- function(lim, nu, rho, sigma2, resol = 1,
                          mean_function = NULL, raster_like = FALSE){
  if (is.null(mean_function)){
    mean_function <- function(z){
      return(0)
    }
  }
  
  xgrid <- seq(lim[1], lim[2], by = resol)
  ygrid <- seq(lim[3], lim[4], by = resol)
  coords <- as.matrix(expand.grid(xgrid, ygrid))
  
  simu <- geoR::grf( grid = coords,  cov.pars = c(sigma2, rho), kappa = nu)
  
  if(!is.null(mean_function)){
    means <- apply(coords, 1, mean_function)
  } else {
    means <- 0
  }
  vals <- means + simu$data
  vals <- vals - min(vals)
  if (raster_like){
    return(list(x = xgrid, y = ygrid, z = matrix(vals, nrow = length(xgrid))))
  }
  return(data.frame(x = coords[, 1], y = coords[, 2], val = vals))
}

#' Transform (x,y,z) list to data.frame
#'
#' @name rasterToDataFrame
#' @param my_list A raster like list containing elements x, y, z
#' (see simSpatialCov)
#' @param level An optional factor argument to add to a level column
#' (useful to pu color for instance)
#' @return A data.frame (convenient for ggplot2)
#' @export
rasterToDataFrame <- function(my_list, level = NULL){
  coords <- as.matrix(expand.grid(my_list$x, my_list$y))
  vals <- as.numeric(my_list$z)
  out <- data.frame(x = coords[, 1], y = coords[, 2], val = vals)
  if (!is.null(level)){
    out$level <- rep(level, nrow(out))
  }
  return(out)
}
