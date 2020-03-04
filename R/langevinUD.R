
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
#' @param alpha Confidence level (default: 0.95, i.e. 95\% confidence intervals)
#' @param leverage Logical. If TRUE, the standardized residuals and the leverage
#' are returned. Default: FALSE. Might not work when there are many observations,
#' because it creates an n times n matrix.
#'
#' @return A list of: est, the vector of estimates, and var, the
#' covariance matrix of the estimates.
#'
#'@export
langevinUD <- function(locs, times, ID = NULL, grad_array, with_speed = TRUE,
                       alpha = 0.95, leverage = FALSE){
  # Check input types
  if (!(inherits(locs, "matrix") & typeof(locs) %in% c("double", "integer")))
    stop("locs must be a numeric matrix")
  if (inherits(grad_array, "matrix")){
    grad_array <- array(grad_array, dim = c(n, 2, 1))
    warning("gradsArray was a matrix, and has been transformed to an array")
  }
  
  # Number of locations
  n <- nrow(locs)
  
  # Check dimensions of inputs
  if (is.null(ID))
    ID <- rep(1, n)
  if (length(times) != n)
    stop("Length of times must be the nrow of locs")
  if (length(ID) != n)
    stop("Length of ID must be the nrow of locs")
  if (dim(grad_array)[1] != n)
    stop("The first dimension of gradientArray must be of size nrow(locs)")
  
  # Number of covariates
  J <- dim(grad_array)[3]
  
  # Find first and last indices of tracks
  break_ind <- which(ID[-1] != ID[-n])
  start_ind <- c(1, break_ind + 1)
  end_ind <- c(break_ind, n)
  
  # Matrix of covariate gradients
  ## doesn't work witha unique covariable
  if(J == 1){
    grad_mat <- 0.5 * matrix(c( grad_array[-end_ind,1,1],  grad_array[-end_ind, 2, 1] ), ncol = 1)
  } else {
    grad_mat <- 0.5 * rbind(grad_array[-end_ind, 1, ], grad_array[-end_ind, 2, ])
  }
  
  # Vector of time differences
  time_lag <- rep(times[-start_ind] - times[-end_ind], 2)
  sq_time_lag <- sqrt(time_lag)
  
  # Derive normalized location increments (response variable)
  loc_increment <- as.numeric(locs[-start_ind, ] - locs[-end_ind, ])
  Y <- loc_increment / sq_time_lag
  
  # (Z'Z)^{-1}
  nu_hat_var <- solve (t(grad_mat * time_lag) %*% grad_mat)
  
  # Degrees of freedom
  DF <- length(Y) - J
  if (DF <= 4){
    stop("nrow(locs) must be strictly larger than 2 + dim(grad_array)[3] / 2")
  }
  
  # Estimate of nu = gamma^2 * beta
  nu_hat <- nu_hat_var %*% t(grad_mat) %*% loc_increment;
  
  # Predicted response
  predictor <- sq_time_lag * grad_mat %*% nu_hat
  
  if (with_speed) {
    # Speed parameter estimate
    gamma2_hat <- colSums( (Y -  predictor) ^ 2 ) / DF
    
    # Habitat selection parameter estimates
    beta_hat <- nu_hat / gamma2_hat * (DF - 2) / DF
    beta_hat_var <- (2 * beta_hat %*% t(beta_hat) / (DF - 4)
                     + nu_hat_var / gamma2_hat * (1 + 2 / (DF - 4)))
  } else {
    # Habitat selection parameter estimates if gamma = 1
    beta_hat <- nu_hat
    beta_hat_var <- nu_hat_var
  }
  
  # Confidence intervals for habitat selection parameters
  quantlow <- (1 - alpha)/2 # Lower quantile
  quantup <- 1 - (1 - alpha)/2 # Upper quantile
  conf_interval_beta <- t(sapply(1:length(beta_hat), function(j) {
    beta_hat[j] + c(1, -1) * stats::qnorm(quantlow) *  sqrt(beta_hat_var[j, j])
  }))
  
  # Confidence interval for speed parameter
  conf_interval_gamma2 <- gamma2_hat * DF / stats::qchisq(c(quantup, quantlow), DF)
  conf_interval <- rbind(conf_interval_beta, conf_interval_gamma2)
  
  # Format output
  rownames(beta_hat_var) <- colnames(beta_hat_var) <- paste0("beta", 1:J)
  rownames(conf_interval) <- c(rownames(beta_hat_var), "gamma2")
  colnames(conf_interval) <- c(quantlow, quantup)
  
  # R squared
  r_square <- 1 - colSums( (Y -  predictor) ^ 2) / sum(Y ^ 2)
  
  # Residuals
  res <- matrix(Y - predictor, ncol = 2)
  
  if(leverage) {
    # Design matrix
    Z <- sq_time_lag * grad_mat
    
    # Hat matrix (leverages on diagonal)
    # (Can be too big when many observations)
    H <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
    
    lever <- diag(H)
    
    # Standardized (Studentized) residuals
    res <- res / (sqrt(gamma2_hat * (1 - diag(H))))
  } else {
    lever <- NULL
  }
  get_predicted <- function(mat_pred){
    purrr::map_dfr(unique(ID), 
                   function(id){
                     my_x <- locs[ID == id, 1]
                     my_y <- locs[ID == id, 2]
                     my_inc_x <- c(mat_pred[ID[-end_ind] == id, 1], NA)
                     my_inc_y <- c(mat_pred[ID[-end_ind] == id, 2], NA)
                     my_times <- times[ID == id]
                     my_sq_time_lag <- c(sqrt(diff(my_times)), NA)
                     tibble::tibble(ID = id, 
                                    x = dplyr::lag(my_x + my_sq_time_lag * my_inc_x),
                                    y = dplyr::lag(my_y + my_sq_time_lag * my_inc_y),
                                    t = my_times)
                   })
  }
  get_increment_residuals <- function(mat){
    purrr::map_dfr(unique(ID), 
                   function(id){
                     ok  <- ID[-end_ind] == id
                     z <- mat[ok, ]
                     loc_times <- times[ID == id]
                     tibble::tibble(ID = id, 
                            x = z[, 1], y = z[, 2]) %>% 
                       dplyr::bind_rows(tibble::tibble(ID = id, x = NA, y = NA)) %>% 
                       dplyr::mutate(t = loc_times)
                   })
  }
  out_res <- get_increment_residuals(res)
  out_pred <- get_predicted(matrix(predictor, ncol = 2))
  
  # Get AIC for fitted model
  AIC <- AICEuler(beta = as.numeric(beta_hat), gamma2 = gamma2_hat,
                  locs = locs, times = times, ID = ID, grad_array = grad_array)
  
  return(list(betaHat = as.numeric(beta_hat), gamma2Hat  = gamma2_hat,
              betaHatVariance = beta_hat_var, CI = conf_interval,
              predicted = out_pred,
              R2 = r_square, residuals = out_res, 
              lever = lever, AIC = AIC))
}
