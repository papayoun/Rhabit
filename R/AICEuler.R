
#' Akaike Information Criterion for Langevin movement model
#' under Euler discretization
#' 
#' @param beta Parameters of the utilisation distribution
#' @param gamma2 Diffusion parameter
#' @param locs Matrix of observed locations (two columns: x, y)
#' @param times Vector of times of observations
#' @param grad_array Three-dimensional array of gradients of covariate fields. 
#' The rows index time, the columns are the dimensions (x and y), and the 
#' layers index the covariates.
#' 
#' @return AIC
AICEuler <- function(beta, gamma2 = 1, locs, times, ID = NULL, grad_array)
{
    # Number of parameters to estimate
    npar <- length(beta)
    if(gamma2 != 1)
        npar <- npar + 1

    # Negative log-likelihood
    nllk <- nllkEuler(beta = beta, gamma2 = gamma2, locs = locs, times = times,
                      ID = ID, grad_array = grad_array)

    AIC <- 2 * (npar + nllk)

    return(AIC)
}
