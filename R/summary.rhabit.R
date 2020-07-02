#' Obtaining the summary of a fitted  Langevin movement model
#'
#' @name summary.rhabit
#' @param object an object of class rhabit, obtained through the function fit_langevin_ud
#' @param ... additional parameters to fit the general usage of the summary function
#'
#' @return the table corresponding to the test of H0 $beta_i = 0 and the table corresponding to the test H0 = "gamma2 =1"
#'
#' @examples
#' data(tracks)
#' fitted_langevin <- fit_langevin_ud(cbind(x,y) ~ grad_c1, data = tracks)
#' summary(fitted_langevin)
#'
#'@export
summary.rhabit <- function(object, ...){
  cat("\nFormula:\n", paste(deparse(object$formula), sep = "\n",
                            collapse = "\n"),
      "\n\n", sep = "")

  ## test for beta =0
  cat("Coefficients:\n")
  est <- object$coefficients
  cov <- object$betaHatVariance
  rdf <- object$df

  se <- sqrt(diag(cov))
  tval <- est/se
  ans_coef <- cbind(Estimate = est, `Std. Error` = se,
               `t value` = tval, `Pr(>|t|)` = 2 * stats::pt(abs(tval), rdf,
                                                     lower.tail = FALSE))
  print(ans_coef)
  ## test for gamma2 =1
  cat("\n\nSpeed parameter:\n")
  gamma2_hat <- object$gamma2
  names(gamma2_hat) <- 'gamma2'
  ans_speed <- cbind(Estimate = gamma2_hat,
                    `chi value` = gamma2_hat*rdf,
                    `Pr(>|t|)` = 2 * min(stats::pchisq(gamma2_hat* rdf, rdf,
                                               lower.tail = TRUE),
                                         stats::pchisq(gamma2_hat*rdf, rdf,
                                    lower.tail = FALSE)) )
  print(ans_speed)
  invisible(list(coeff = ans_coef, speed = ans_speed))
}


#' Show the estimated speed of the Langevin model
#'
#' @param x a rhabit object obtained by fit_langevin_ud
#'
#' @return the diffusion parameter
#' @export
#'
#' @examples
#' data(tracks)
#' fitted_langevin <- fit_langevin_ud(cbind(x,y) ~ grad_c1, data = tracks)
#' speed_coef(fitted_langevin)
speed_coef <- function(x){
  x$gamma2
}

#' Extract the coefficients associated with the covariates in a rhabit model
#'
#' @param object a rhabit model fit using fit_langevin_ud
#' @param ... additional arguments. 
#' @param complete argument specifies whether  NA should be displayed or not. its value is set to TRUE by default
#' 
#' @return the coefficients of the fitted rhabit  model
#' @export
#'
#' @examples
#' data(tracks)
#' fitted_langevin <- fit_langevin_ud(cbind(x,y) ~ grad_c1, data = tracks)
#' coef(fitted_langevin)

coef.rhabit <- function(object, complete = TRUE, ...){
  cf <- object$coefficients
  if (complete)
    cf
  else cf[!is.na(cf)]
}



#' Print the formula used to specify the model and the estimated coefficients
#'
#' @param x  a rhabit model fit using fit_langevin_ud
#' @param ... additionnal arguments, might be digits a non-null value specifying the minimum number 
#' of digits to be printed in values.
#'
#' @return The rhabit object passed in argument is returned.
#' @export
#'
#' @examples
#' data(tracks)
#' fitted_langevin <- fit_langevin_ud(cbind(x,y) ~ grad_c1, data = tracks)
#' print(fitted_langevin)

print.rhabit <- function(x, ...){
  if (! 'digits' %in% names(list(...)) )
    digits <- 3
  cat("\nFormula:\n", paste(deparse(x$formula), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(coef.rhabit(x))) {
    cat("Coefficients:\n")
    print.default(format(coef.rhabit(x), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No coefficients\n")

  if (length(speed_coef(x))) {
    cat("Speed:\n")
    print.default(format(speed_coef(x), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }

    cat("\n")
  invisible(x)
  cat("\n ")

}
