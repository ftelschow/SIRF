#------------------------------------------------------------------------------#
#                                                                              #
#     Functions to compute parametric quantiles
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#
#
#------------------------------------------------------------------------------#
# Developer notes:
#
#------------------------------------------------------------------------------#
#' This functions computes the quatiles of the maximum of combinations of
#' statistics of the form -X[i], X[i] or |X[i]|  where for {1,...,I} the
#' random variables X[i] are independently standard Gaussian distribtued.
#'
#' @inheritParams PreimageC
#' @param p numeric probability for which the quantile should be computed.
#' @param muC list output from PreimageC
#' @return quantile of the maximum of the random variables
#' @export
maxGauss_quantile <- function(p, muC){
  Mpm   = sum(muC$minus && muC$plus)
  M     = sum(muC$minus) + sum(muC$plus) - 2 * Mpm

  if(Mpm == 0 && M != 0){
    q <- qnorm(p^(1 / M))
  }else if(Mpm != 0 && Mpm == M){
    q <- VGAM::qfoldnorm(p^(1 / Mpm ))
  }else if(M != 0){
    # Find the quantile
    minFunction <- function(q){
      # remove dependency on VGAM
      pnorm(q) * VGAM::pfoldnorm(q)^(Mpm/M) - p^(1/M)
    }
    q = uniroot(minFunction, interval = c(-50,50), tol = 1e-5)$root
  }else{
    q = 0
  }
  return(q)
}

#' This functions computes the probability to be below q for the maximum of
#' combinations of statistics of the form -X[i], X[i] or |X[i]|  where for
#' {1,...,I} the random variables X[i] are independently standard Gaussian
#' distribtued.
#'
#' @inheritParams PreimageC
#' @param p numeric probability for which the quantile should be computed.
#' @param muC list output from PreimageC
#' @return quantile of the maximum of the random variables
#' @export
maxGauss_p <- function(q, muC){
  Mpm   = sum(muC$minus && muC$plus)
  M     = sum(muC$minus) + sum(muC$plus) - 2 * Mpm

  if(Mpm == 0 && M != 0){
    p <- pnorm(q)^M
  }else if(Mpm != 0 && Mpm == M){
    p <- VGAM::pfoldnorm(q)^Mpm
  }else if(M != 0){
    p = pnorm(q)^M * VGAM::pfoldnorm(q)^Mpm
  }else{
    p = 1
  }
  return(p)
}

#' This functions computes the quatiles of the maximum of combinations of
#' statistics of the form -X[i], X[i] or |X[i]|  where for {1,...,I} the
#' random variables X[i] are independently Student's t-distribtued with df
#' degrees of freedom
#'
#' @inheritParams PreimageC
#' @param p numeric probability for which the quantile should be computed.
#' @param muC list output from PreimageC
#' @param df numeric an integer giving the degrees of freedom
#' @return quantile of the maximum of the random variables
#' @export
maxT_quantile <- function(p, muC, df){
  Mpm   = sum(muC$minus & muC$plus)
  M     = sum(muC$minus) + sum(muC$plus) - Mpm

  if(Mpm == 0 && M != 0){
    q <- qt(p = p^(1 / M), df = df)
  }else if(Mpm != 0 && Mpm == M){
    q <- extraDistr::qht(p = p^(1 / Mpm ), nu = df)
  }else if(M != 0){
    # Find the quantile by finding a root
    minFunction <- function(q){
      # remove dependency on VGAM
      pt(q = q, df = df) * extraDistr::pht(q = q, nu = df)^(Mpm/M) - p^(1/M)
    }
    q = uniroot(minFunction, interval = c(-50,50), tol = 1e-6)$root
  }else{
    q = 0
  }
  return(q)
}

#' This functions computes the quatiles of the maximum of combinations of
#' statistics of the form -X[i], X[i] or |X[i]|  where for {1,...,I} the
#' random variables X[i] are independently Student's t-distribtued with df
#' degrees of freedom
#'
#' @inheritParams PreimageC
#' @param p numeric probability for which the quantile should be computed.
#' @param muC list output from PreimageC
#' @param df numeric an integer giving the degrees of freedom
#' @return quantile of the maximum of the random variables
#' @export
maxT_p <- function(q, muC, df){
  Mpm   = sum(muC$minus & muC$plus)
  M     = sum(muC$minus) + sum(muC$plus) - Mpm

  if(Mpm == 0 && M != 0){
    p <- pt(q = q, df = df)^M
  }else if(Mpm != 0 && Mpm == M){
    p <- extraDistr::pht(q = q, nu = df)^Mpm
  }else if(M != 0){
    # Find the quantile by finding a root
    minFunction <- function(q){
      # remove dependency on VGAM
      pt(q = q, df = df)^M * extraDistr::pht(q = q, nu = df)^Mpm
    }
    p = uniroot(minFunction, interval = c(-50,50), tol = 1e-6)$root
  }else{
    p = 0
  }
  return(p)
}
