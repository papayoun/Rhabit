
#' Rhabit to raster
#' 
#' Transform a (x,y,z) list, as used in Rhabit, into a raster object
#' 
#' @param r List of three elements: x (grid of x values), y (grid of y values),
#' and z (matrix of values)
#' 
#' @return Raster object
#' 
#' @export
rhabitToRaster <- function(r) {
  raster::rasterFromXYZ(rasterToDataFrame(r))
}
