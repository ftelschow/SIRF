% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoPEFunctions.R
\name{sim_SCoPES}
\alias{sim_SCoPES}
\title{This functions computes the SCoPES corresponding to an estimator and a set
of functions given as a matrix with columns being the cut-off functions.}
\usage{
sim_SCoPES(
  Msim,
  N,
  alpha,
  C,
  q_method,
  model,
  I = NULL,
  inclusion = list(L = "inequal", U = "inequal")
)
}
\arguments{
\item{alpha}{numeric between 0 and 1. It will produce (1-alpha) SCoPES.}

\item{C}{T x nlvl matrix containing nlvl functions for which the preimage
should be computed. At the moment only the first and last function
matter, i.e., they define a tube.}

\item{inclusion}{list with elements
 \itemize{
  \item L string either "inequal" or "equal" determining whether the lower
  excursion set uses < c or <=c.
  \item U string either "inequal" or "equal" determining whether the upper
  excursion set uses > c or >=c.
}}
}
\value{
Standard error under the assumption the data is Gaussian
}
\description{
This functions computes the SCoPES corresponding to an estimator and a set
of functions given as a matrix with columns being the cut-off functions.
}
