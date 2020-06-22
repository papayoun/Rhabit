#' Obtaining the summaryof a fitted  Langevin movement model 
#'
#' @name summary.rhabit
#' @param rhabit_obj an object of class rhabit, obtained through the function langevinUD
#' 
#'
#' @return the table corresponding to the test of H0 $beta_i = 0$ 
#'
#'@export 
summary.rhabit <- function(rhabit_obj){
  est <- rhabit_obj$coefficients
  cov <- rhabit_obj$betaHatVariance
  rdf <- rhabit_obj$df
  
  se <- sqrt(diag(cov))
  tval <- est/se
  ans <- cbind(Estimate = est, `Std. Error` = se, 
               `t value` = tval, `Pr(>|t|)` = 2 * pt(abs(tval), rdf, 
                                                     lower.tail = FALSE))
  return(ans)
} 



coef.rhabit <- function(object, complete = TRUE){
  cf <- object$coefficients
  if (complete) 
    cf
  else cf[!is.na(cf)]
}



#' Title
#'
#' @param rhabit_obj 
#'
#' @return
#' @export
#'
#' @examples
print.rhabit <- function(x,digits = max(3L, getOption("digits") - 3L)){
  
  cat("\nCall:\n", paste(deparse(x$formula), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L, 
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
  cat("\n ")
  
} 
