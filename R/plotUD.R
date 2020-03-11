#' Plotting the classical RSF UD
#' @name plotUD
#'
#' @param covariates covariates list, that must be a list of "raster-like"
#' elements. It is for the moment supposed that the grid is regular and
#' shared by all elemebts of covariates
#' @param beta beta values, must be the same length as covariates
#' @param log logical: get log value? by default = F
#'
#' @return a ggplot
#' @export
plotUD <- function(covariates, beta, log = F){
  ud_df <- getUD(covariates, beta) %>% 
    rasterToDataFrame() 
  ggplot(ud_df, 
         aes(x = x, y = y)) +
    geom_raster(aes(fill = val)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_viridis_c() +
    labs(fill = expression(pi(z)))
}