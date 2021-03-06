% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/val.R
\name{val}
\alias{val}
\title{Estimate the coefficients}
\usage{
val(x, z, y, p, nknots, m = 2, lambda = NULL, cv = FALSE,
  sigma_error = diag(1, length(y)))
}
\arguments{
\item{x}{predictors, a vector}

\item{z}{assist variable, a vector}

\item{y}{dependent variable, a vector}

\item{p}{order of B-splines}

\item{m}{order of penalty term, m = 2 by default}

\item{lambda}{coefficient of penalty term, need to be specified}

\item{sigma_error}{covariance matrix of disturbance term u}

\item{nkonts}{the number of knots}
}
\description{
Estimate the coefficients
}

