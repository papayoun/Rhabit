#' Obtaining the classical RSF UD
#' @name getUD
#'
#' @param covariates covariates list, that must be a list of "raster-like"
#' elements. It is for the moment supposed that the grid is regular and
#' shared by all elemebts of covariates
#' @param beta beta values, must be the same length as covariates
#' @param log logical: get log value? by default = F
#'
#' @return a raster-like list
#' @export
getUD <- function(covariates, beta, log = F){
  ud_rast <- covariates[[1]] # initialization
  dx <- diff(ud_rast$x)[1]
  dy <- diff(ud_rast$y)[1]
  J <- length(covariates)
  ud_rast$z <- Reduce("+",
                      lapply(1:J, function(j)
                        dx * dy * beta[j] * covariates[[j]]$z))
  if (!log){
    ud_rast$z <- exp(ud_rast$z)
    ud_rast$z <- ud_rast$z / sum(ud_rast$z)# Normalization
  }
  ud_rast
}
