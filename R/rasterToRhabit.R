
#' Raster to Rhabit
#' 
#' Transform a raster into a (x,y,z) list, as used in Rhabit
#' 
#' @param r Raster object
#' 
#' @return List of three elements: x (grid of x values), y (grid of y values),
#' and z (matrix of values)
#' 
#' @importFrom raster extent res as.matrix
#' 
#' @export
rasterToRhabit <- function(r) {
  # Extract limits and resolution from raster
  lim <- as.vector(extent(r))
  res <- res(r)
  # Define x and y grids
  xgrid <- seq(lim[1] + res[1]/2, lim[2] - res[1]/2, by=res[1])
  ygrid <- seq(lim[3] + res[2]/2, lim[4] - res[2]/2, by=res[2])
  
  # Put values in the right matrix format
  z <- t(apply(as.matrix(r),2,rev))
  
  return(list(x = xgrid, y = ygrid, z = z))
}
