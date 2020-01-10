## find neighbours around x0

getGridZoom <- function(covar, x0, lag_inter = 2){
  if (is.null(covar))
    return(NULL)
  else{
    x_pos <- findInterval(x0[1], covar$x)
    y_pos <- findInterval(x0[2], covar$y)
    x_inter <- max(1, x_pos - lag_inter):min(length(covar$x), x_pos + lag_inter)
    y_inter <- max(1, y_pos - lag_inter):min(length(covar$y), y_pos + lag_inter)
    return(list(x = covar$x[x_inter],
                y = covar$y[y_inter],
                z = covar$z[x_inter, y_inter]))
  }
}
