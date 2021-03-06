% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.rhabit.R
\name{summary.rhabit}
\alias{summary.rhabit}
\title{Obtaining the summary of a fitted  Langevin movement model}
\usage{
\method{summary}{rhabit}(object, ...)
}
\arguments{
\item{object}{an object of class rhabit, obtained through the function fit_langevin_ud}

\item{...}{additional parameters to fit the general usage of the summary function}
}
\value{
the table corresponding to the test of H0 $beta_i = 0 and the table corresponding to the test H0 = "gamma2 =1"
}
\description{
Obtaining the summary of a fitted  Langevin movement model
}
\examples{
data(tracks)
fitted_langevin <- fit_langevin_ud(cbind(x,y) ~ grad_c1, data = tracks)
summary(fitted_langevin)

}
