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
  ggplot2::ggplot(ud_df, 
                  ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$val)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(fill = expression(pi(z)))
}