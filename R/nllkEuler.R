
#' Negative log-likelihood function for Langevin movement model
#' under Euler discretization
#' 
#' @param beta Parameters of the utilisation distribution
#' @param gamma2 Diffusion parameter
#' @param locs Matrix of observed locations (two columns: x, y)
#' @param times Vector of times of observations
#' @param ID Vector of track identifiers, must be of length n
#' @param grad_list List of gradients of covariates, evaluated at
#' the locations of the observations, must be 3d array of dim(n, 2, J)
#' where J is the amount of covariates .
#' @keywords internal
#' @return Negative log-likelihood
nllkEuler <- function(beta, gamma2 = 1, locs, times, ID = NULL, grad_list) {
    n <- nrow(locs)
    # multiply gradients by beta coefficients
    nu <- .5 * beta * gamma2 
    gradmat_list <- lapply(seq_len(length(grad_list)),  function(d_){
        grad_list[[d_]] * nu[d_]
    })
    gradmat <- Reduce('+', gradmat_list)
     
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
