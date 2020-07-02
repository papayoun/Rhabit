#' Obtaining the summary of a fitted  Langevin movement model
#'
#' @name summary.rhabit
#' @param object an object of class rhabit, obtained through the function fit_langevin_ud
#'
#'
#' @return the table corresponding to the test of H0 $beta_i = 0 and the table corresponding to the test H0 = "gamma2 =1"
#'
#' @examples
#' data(tracks)
#' fitted_langevin <- fit_langevin_ud(cbind(x,y) ~ grad_c1, data = tracks)
#' summary(fitted_langevin)
#'
#'@export
summary.rhabit <- function(object){
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
               `t value` = tval, `Pr(>|t|)` = 2 * pt(abs(tval), rdf,
                                                     lower.tail = FALSE))
  print(ans_coef)
  ## test for gamma2 =1
  cat("\n\nSpeed parameter:\n")
  gamma2_hat <- object$gamma2
  names(gamma2_hat) <- 'gamma2'
  ans_speed <- cbind(Estimate = gamma2_hat,
                    `chi value` = gamma2_hat*rdf,
                    `Pr(>|t|)` = 2 * min(pchisq(gamma2_hat* rdf, rdf,
                                               lower.tail = TRUE),
                                        pchisq(gamma2_hat*rdf, rdf,
                                    lower.tail = FALSE)) )
  print(ans_speed)
  invisible(list(coeff = ans_coef, speed = ans_speed))
}


#' Show the estimated speed of the Langevin model
#'
#' @param x a rhabit object obtained by fit_langevin_ud
#'
#' @return
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
#' @param complete  should NA be shown
#'
#' @return the coefficients of the fitted rhabit  model
#' @export
#'
#' @examples
#' data(tracks)
#' fitted_langevin <- fit_langevin_ud(cbind(x,y) ~ grad_c1, data = tracks)
#' coef(fitted_langevin)

coef.rhabit <- function(object, complete = TRUE){
  cf <- object$coefficients
  if (complete)
    cf
  else cf[!is.na(cf)]
}



#' Title
#'
#' @param x  a rhabit model fit using fit_langevin_ud
#' @param digits a non-null value for digits specifies the minimum number o
#' f significant digits to be printed in values.
#'
#' @return
#' @export
#'
#' @examples
#' data(tracks)
#' fitted_langevin <- fit_langevin_ud(cbind(x,y) ~ grad_c1, data = tracks)
#' print(fitted_langevin)

print.rhabit <- function(x,digits = max(3L, getOption("digits") - 3L)){
  cat("\nFormula:\n", paste(deparse(x$formula), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L,
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
