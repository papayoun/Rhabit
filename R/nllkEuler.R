
#' Negative log-likelihood function for Langevin movement model
#' under Euler discretization
#' 
#' @param beta Parameters of the utilisation distribution
#' @param gamma2 Diffusion parameter
#' @param locs Matrix of observed locations (two columns: x, y)
#' @param times Vector of times of observations
#' @param ID Vector of track identifiers, must be of length n
#' @param grad_array Three-dimensional array of gradients of covariate fields. 
#' The rows index time, the columns are the dimensions (x and y), and the 
#' layers index the covariates.
#' 
#' @return Negative log-likelihood
nllkEuler <- function(beta, gamma2 = 1, locs, times, ID = NULL, grad_array)
{
    n <- nrow(locs)
    
    # multiply gradients by beta coefficients
    gradmat <- 0.5 * gamma2 * apply(grad_array, 2, function(mat) mat %*% beta)
    
    # only one track
    if(is.null(ID))
        ID <- rep(1,nrow(locs))
    
    # indices of first and last obs in tracks
    break_ind <- which(ID[-1] != ID[-n])
    start_ind <- c(1, break_ind + 1) # first obs
    end_ind <- c(break_ind, n) # last obs

    # time intervals
    dt <- times[-start_ind] - times[-end_ind]

    # log likelihood
    llk <- sum(dnorm(locs[-start_ind,],
                     locs[-end_ind,] + dt * gradmat[-end_ind,],
                     sqrt(gamma2 * dt),
                     log=TRUE))
    
    return(-llk)
}
