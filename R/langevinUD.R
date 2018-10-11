#' Obtaining UD estimate of an animal using Langevin movement model
#'
#' @name langevinUD
#' @param locs Matrix of locations having size, say, n
#' @param times Vector of observation times, must be of length n
#' @param ID Vector of track identifiers, must be of length n
#' @param grad_array Array of gradients of covariates, evaluated at
#' the locations of the observations, must be 3d array of dim(n, 2, J)
#' where J is the amount of covariates
#' @param with_speed Logical. If TRUE, the speed parameter is estimated
#' Other wise it is set to one
#'
#' @return A list of: est, the vector of estimates, and var, the
#' covariance matrix of the estimates.
#'
#'
langevinUD <- function(locs, times, ID = NULL, grad_array, with_speed = TRUE){
  if (!(inherits(locs, "matrix") & typeof(locs) %in% c("double", "integer")))
    stop("locs must be a numeric matrix")
  n <- nrow(locs)
  if (inherits(grad_array, matrix)){
    grad_array <- array(grad_array, dim = c(n, 2, 1))
    warning("gradsArray was a matrix, and has been transformed to an array")
  }
  J <- dim(grad_array)[3]
  if (length(times) != n)
    stop("Length of times must be the nrow of locs")
  if (is.null(ID))
    ID <- rep(1, n)
  if (length(ID) != n)
    stop("Length of ID must be the nrow of locs")
  if (dim(grad_array)[1] != n)
    stop("The first dimension of gradientArray must be of size nrow(locs)")
  break_ind <- which(ID[-1] != ID[-n])
  start_ind <- c(1, break_ind + 1)
  end_ind <- c(break_ind, n)
  grad_mat <- 0.5 * rbind(grad_array[-end_ind, 1, ], grad_array[-end_ind, 2, ])
  time_lag <- rep(times[-start_ind] - times[-end_ind], 2)
  sq_time_lag <- sqrt(time_lag)
  loc_increment <- as.numeric(locs[-start_ind, ] - locs[-end_ind, ])
  Y <- loc_increment / sq_time_lag
  nu_hat_var <- solve (t(grad_mat * time_lag) %*% grad_mat)
  DF <- length(Y) - J
  if(DF <= 4){
    stop("nrow(locs) must be strictly larger than 2 + dim(grad_array)[3] / 2")
  }
  nu_hat <- nu_hat_var %*% t(grad_mat) %*% loc_increment;
  gamma2_hat <- NULL
  if (with_speed){
    predictor <- sq_time_lag * grad_mat %*% nu_hat
    gamma2_hat <- colSums( (Y -  predictor) ^ 2 ) / DF
    beta_hat <- nu_hat / gamma2_hat * (DF - 2) / DF
    beta_hat_var <- (2 * beta_hat %*% t(beta_hat) / (DF - 4)
                        + nu_hat_var / gamma2_hat * (1 + 2 / (DF - 4)))
  }
  else{
    beta_hat <- nu_hat
    beta_hat_var <- nu_hat_var
  }
  conf_interval <- t(sapply(1:length(beta_hat), function(j){
    beta_hat[j] + c(1, -1) * stats::qnorm(0.025) *  sqrt(beta_hat_var[j, j])
  }))
  rownames(beta_hat_var) <- colnames(beta_hat_var) <- paste0("beta", 1:J)
  rownames(conf_interval) <- rownames(beta_hat_var)
  return(list(betaHat = as.numeric(beta_hat), gamma2Hat  = gamma2_hat,
              betaHatVariance = beta_hat_var, betaHat95CI = conf_interval))
}
