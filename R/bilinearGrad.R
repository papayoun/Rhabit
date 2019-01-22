
#' Extract covariate grid cell
#' 
#' @param loc Location
#' @param xgrid Grid along x axis
#' @param ygrid Grid along y axis
#' @param covmat Covariate values on (xgrid,ygrid)
#' 
#' @return List of two elements: coords (coordinates of bottom-left and
#' top-right corners of the grid cell), and values (values of the covariate 
#' at the four corners of the grid cell).
gridCell <- function(loc, xgrid, ygrid, covmat)
{
  ix <- findInterval(loc[1], xgrid)
  iy <- findInterval(loc[2], ygrid)
  
  # coordinates of bottom-left and top-right points of grid cell
  coords <- c(xgrid[ix], xgrid[ix+1], ygrid[iy], ygrid[iy+1])
  
  # values of the covariate at the four corner points
  values <- covmat[ix:(ix+1), iy:(iy+1)]
  
  return(list(coords=coords, values=values))
}

#' Gradient for bilinear interpolation
#' 
#' @param loc Point at which the gradient should be evaluated (2d vector)
#' @param cov_list List of J (number of covariates) "raster like" elements.
#' A raster like element is a 3 elements list with named elements
#' 1) "x" a vector of increasing x locations (at which the covariate is sampled)
#' 2) "y" a vector of increasing y locations (at which the covariate is sampled)
#' 3) "z" a size(x)*size(y) matrix giving covariate values at location (x, y)
#' 
#' @return Gradient of all covariates at the given location
#' 
#' @export
bilinearGrad <- function(loc, cov_list) {
  J <- length(cov_list)
  grad_val <- sapply(1:J, function(j) {
    # Grid over x axis
    x_grid <- cov_list[[j]]$x
    # Grid over y axis
    y_grid <- cov_list[[j]]$y
    # Matrix of covariate values
    cov_mat <- cov_list[[j]]$z
    
    # Extract values of coordinate around location of interest
    cell <- gridCell(loc=loc, xgrid=x_grid, ygrid=y_grid, covmat=cov_mat)
    x <- cell$coords[1:2]
    y <- cell$coords[3:4]
    f <- cell$values

    # Gradient along x axis
    dfdx <- 
      ((y[2] - loc[2]) * (f[2,1] - f[1,1]) + (loc[2] - y[1]) * (f[2,2] - f[1,2])) /
      ((y[2] - y[1]) * (x[2] - x[1]))
    
    # Gradient along y axis
    dfdy <- 
      ((x[2] - loc[1]) * (f[1,2] - f[1,1]) + (loc[1] - x[1]) * (f[2,2] - f[2,1])) /
      ((y[2] - y[1]) * (x[2] - x[1]))
    
    return(c(dfdx, dfdy))
  })
  
  return(grad_val)
}

#' Apply bilinearGrad to matrix of locations
#' 
#' @param locs Matrix of locations (with two columns)
#' @param cov_list List of J (number of covariates) "raster like" elements.
#' A raster like element is a 3 elements list with named elements
#' 1) "x" a vector of increasing x locations (at which the covariate is sampled)
#' 2) "y" a vector of increasing y locations (at which the covariate is sampled)
#' 3) "z" a size(x)*size(y) matrix giving covariate values at location (x, y)
#' 
#' @return An array of dimension (n, 2, J), where n is the number of locations,
#' and J is the number of covariates.
#' 
#' @export
bilinearGradArray <- function(locs, cov_list) {
  # Compute gradient at each location
  grad <- unlist(lapply(1:nrow(locs), function(i) 
    bilinearGrad(loc = locs[i,], cov_list = cov_list)))
  
  # Put into array with correct dimensions
  gradarray <- array(grad, c(2, length(cov_list), nrow(locs)))
  
  # Swap dimensions to obtain an array with the correct format
  gradarray <- aperm(gradarray, c(3, 1, 2))
  
  return(gradarray)
}
