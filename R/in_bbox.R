#' Check if a point is in the bouding box
#'
#' @param point The point to be checked a c(x,y) vector
#' @param bbox The bounding box bbox[i,j] if i=1 is the minimal value for the bbox and i=2 the maximal value, j=1 corresponds to the x dimension and j=2 to the y dimension  
#' @keywords internal
#' @return a boolean
 in_bbox <- function(point, bbox){
   return( bbox[1,1] <= point[1] && bbox[2,1] >= point[1]  && bbox[1,2] <= point[2] && bbox[2,2] >= point[2] ) 
 }