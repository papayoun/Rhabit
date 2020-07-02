
#' Akaike Information Criterion for Langevin movement model
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
#' 
#' 
#' @return the approximated AIC using Euler discretization scheme
#' 
#' @export
AICEuler <- function(beta, gamma2 = 1, locs, times, ID = NULL, grad_list)
{
    # Number of parameters to estimate
    npar <- length(beta)
    if(gamma2 != 1)
        npar <- npar + 1
    # Negative log-likelihood
    nllk <- nllkEuler(beta = beta, gamma2 = gamma2, locs = locs, times = times,
                      ID = ID, grad_list = grad_list)

    AIC <- 2 * (npar + nllk)

    return(AIC)
}
