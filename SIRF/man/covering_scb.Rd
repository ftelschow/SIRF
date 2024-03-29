% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ConfidenceBands.R
\name{covering_scb}
\alias{covering_scb}
\title{Computes the quantile of the maximum statistic}
\usage{
covering_scb(
  Msim = 50000,
  Nvec,
  level,
  method = list(name = "GKF"),
  x,
  mu,
  sigma,
  noise,
  transformation,
  moments = NULL,
  bias.est,
  se.est,
  trueValue,
  smoothing = NULL,
  ...
)
}
\arguments{
\item{Msim}{number of simulations.}

\item{Nvec}{vector containing the sample sizes}

\item{level}{numeric the targeted covering probability. Default 0.95.}

\item{method}{string specifying the method to construct the scb, i.e.
estimatate the quantile.
Current options are "tGKF", "GKF", "NonparametricBootstrap", "MultiplierBootstrap".
Default value is "tGKF".}

\item{transformation}{expression of the form f ~ f(mu1, ..., muk ) giving
the transformation of the moments. Can be a "linear", "cohensd", "skewness" or
"kurtosis". Then standard estimators are used.
Current options are "skewness" and "kurtosis".}

\item{moments}{string of the form c("mu1", ..., "muk") containing the names of
the moments contained in the order of transformation input.}

\item{bias.est}{a string indicating the bias estimation method.
 \itemize{
  \item "asymptotic gaussian" means bias = 0
  \item "exact gaussian" uses the true finite sample variance bias assuming
the transdormation is cohensd, skewness or kurtosis
  \item "estimate" means that the bias is estimated from the delta residuals
}
Default is bias = 0, i.e., "asymptotic gaussian".
,}

\item{se.est}{\itemize{
  \item "asymptotic gaussian" means true Gaussian asymptotic variance is used in the
  standard error, i.e. se("skewness") = 6 / sqrt(sample size)
  \item "exact gaussian" uses the true finite sample variance bias assuming
the transdormation is cohensd, skewness or kurtosis
  \item "estimate" means that the bias is estimated from the delta residuals
}}

\item{smoothing}{either NULL indicating "no smoothing applied" or a list.
In the latter case locpoly from KernSmooth is used
 \itemize{
  \item If a field "bandwidth" is contained the value is used for the
  smoothing. Othererwise crossvalidation estimates the optimal bandwidth
  \item  If a field "degree" is contained the value is used for the
  smoothing. Othererwise degree = 1 is used.
   }}
}
\value{
quantile of the maximum process
}
\description{
Computes the quantile of the maximum statistic
}
