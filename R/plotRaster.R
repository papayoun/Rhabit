
#' Plot raster with ggplot
#' 
#' @param rast Raster object
#' @param norm Logical. Should the raster be normalized to integrate to 1?
#' @param log Logical. Should the raster be plotted on the log scale?
#' @param scale.name Name of the scale
#' @param light Logical. If TRUE, the plot is minimalistic.
#' 
#' @importFrom ggplot2 ggplot geom_raster coord_equal geom_point geom_path aes_string
#' theme element_blank
#' @importFrom viridis scale_fill_viridis
#' @importFrom raster xres yres values 
#' @importFrom sp coordinates
#' 
#' @export
plotRaster <- function(rast, norm=FALSE, log=FALSE, scale.name="", light=FALSE)
{
  # Define data frame from raster object
  covmap <- data.frame(coordinates(rast),val=values(rast))
  
  # Normalize and/or take log scale if necessary
  if(norm) {
    s <- sum(covmap$val) * xres(rast) * yres(rast)
    covmap$val <- covmap$val/s
  }
  if(log) {
    covmap$val <- log(covmap$val)
  }
  
  # Creat plot to return
  p <- ggplot(covmap, aes_string(x="x",y="y")) + geom_raster(aes_string(fill="val")) +
    coord_equal() 
  
  # If light, remove scale guide, and axis labels
  if(light) {
    p <- p + scale_fill_viridis(guide = "none") +
      theme(axis.title = element_blank(), axis.text = element_blank(), 
            axis.ticks = element_blank())
  } else {
    p <- p + scale_fill_viridis(name = scale.name)
  }
  
  return(p)
}
