% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R, R/with.qibm.R
\docType{methods}
\name{LRM}
\alias{LRM}
\alias{LRM,qibm-method}
\alias{describe,qibmLRM-method}
\title{Logistic Regression Model Estimates}
\usage{
LRM(x, ...)

\S4method{LRM}{qibm}(x, coef, N, seed = 123, ...)

\S4method{describe}{qibmLRM}(x, alpha = 0.05, alternative = c("two.sided",
  "less", "greater"), scale = 1)
}
\arguments{
\item{x}{an \R object.}

\item{...}{further arguments passed to \code{\link{with}}.}

\item{coef}{two-element vector specifying the logistic regression
intercept and slope with which to simulate outcomes.}

\item{N}{scalar or vector of simulated sample sizes.}

\item{seed}{random number seed for the simulations.}

\item{alpha}{significance for statistical testing of the odds ratio.}

\item{alternative}{direction of the alternative hypothesis.}

\item{scale}{value at which to compute the odds ratio}
}
\description{
Logistic Regression Model Estimates
}
\seealso{
\code{\link{LRMCoef}}.
}
