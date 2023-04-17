#------------------------------------------------------------------------------#
#                                                                              #
#     Insignificance values
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#
#
#------------------------------------------------------------------------------#
# Developer notes:
#
#------------------------------------------------------------------------------#
#' This functions computes some possible insignificance values for the simple iid
#' model.
#'
#' @inheritParams PreimageC
#' @param SCoPEStype vector of same length as mu giving the possible alternative
#' @param Delta
#' @return quantile of the maximum of the random variables
#' @export
IVobs <- function(scopes, method){
  #
  dC = dim(scopes$C)
  if(dC[2]==1){
    scopes$C = cbind(scopes$C,scopes$C)
  }else if(dC[2]>2){
    stop("Currently only one tube is supported.")
  }

  if(method$name == "gauss"){
    return(IVs_IID(scopes))
  }else if(method$name == "t"){
    return(IVs_IID(scopes,
                   cdf = pt,
                   cdf_abs = function(q) extraDistr::pht(q, nu = method$df),
                   df = method$df))
  }else if(method$name == "mboot"){
    # # Get a bootstrap sample
    # mBoot <- MultiplierBootstrapSplit(R = method$R,
    #                                   rep(TRUE, length(scopes$hatmu)),
    #                                   rep(TRUE, length(scopes$hatmu)),
    #                                   alpha   = 0.05,
    #                                   Mboots  = method$Mboots,
    #                                   method  = method$Boottype,
    #                                   weights = method$weights)
    # # Compute the local insignificance values
    # IV1 <- rep(NaN, length(scopes$hatmu))
    # if(sum(index_) != 0){
    #   IV1[scopes$hatLC[, 1]] <- apply(C[scopes$hatLC[, 1], 1] - mBoot$samples[scopes$hatLC[, 1],] <=
    #                                     scopes$hatmu[scopes$hatLC[, 1]], 1, mean )
    # }
    # if(sum(scopes$hatUC[, 2]) != 0){
    #   IV1[scopes$hatUC[, 2]] <- apply(C[scopes$hatUC[, 2], 2] + mBoot$samples[scopes$hatUC[,2],] >=
    #                                     scopes$hatmu[scopes$hatUC[,2]], 1, mean )
    # }
    #
    # # Compute the extreme values
    # IV2 <- 1 - mean(apply(abs(mBoot$samples[index,]) < abs(hatmu[index]), 2, all))
  }else{
    stop("Currently not implemented for the method you request.")
  }
}


#' This functions computes the observational insignificance value (p-value )
#' for an IID Gaussian model the maximum of combinations of statistics of the form
#' -X[i], X[i] or |X[i]| where for {1,...,I} the random variables X[i] are
#' independently standard Gaussian distribtued.
#'
#' @inheritParams PreimageC
#' @param SCoPEStype vector of same length as mu giving the possible alternative
#' @param cdf a cummulative distribution function taking a quantile q as input and
#' computing the probability that the random variable belonging to the cdf stays
#' below q
#' @param ... further input variables for cdf
#' @return quantile of the maximum of the random variables
#' @export
IVs_IID <- function(scopes, cdf = pnorm, cdf_abs = VGAM::pfoldnorm, rcdf = rnorm, Msim = 5e5, KK = 5,...){
  # Global probability of detecting at least one element
  IV1 <- 1 - cdf_abs(scopes$q)^length(scopes$x)

  # The scaled distance process
  if(dim(scopes$C)[2] == 1){
    C = cbind(scopes$C, scopes$C)
  }else{
    C = scopes$C
  }
  TC <- (scopes$hatmu - C) / (scopes$tN * scopes$hatsigma)

  pvals = 1 - cdf_abs(abs(TC[,1]))
  hatm0 = min(2*sum(pvals >= 0.5), length(pvals))

  detectN = sum(scopes$hatUC | scopes$hatLC)
  Tnull <- abs(scopes$hatmu / (scopes$tN * scopes$hatsigma))
  Tnull = sort(Tnull[scopes$hatUC | scopes$hatLC])
  Tnull = Tnull[1:KK]
  # Get the "p-values" for each location s
  IVloc <- rep(NA, length(scopes$hatmu))
  names(IVloc)  <- "C"
  names(IVloc)[scopes$hatLC] <- "L"
  names(IVloc)[scopes$hatUC] <- "U"

  pL <- 1 - cdf(-TC[, 1], ...)
  pU <- 1 - cdf( TC[, 2], ...)

  IVloc <- vapply(1:length(scopes$x), function(x) min(pL[x], pU[x]), FUN.VALUE = 0.1)

  # How extreme are the detections given the null
  IVobs <- prod(ifelse(any(scopes$hatLC), IVloc[scopes$hatLC], 1)) *
                  prod(ifelse(any(scopes$hatUC), IVloc[scopes$hatUC], 1))

  rsample = abs(matrix(rcdf(length(scopes$x) * Msim, ...), nrow = length(scopes$x)))
  IVk <- mean(apply(rsample, 2, function(col) sum(col > scopes$q) >= sum(scopes$hatUC) ) &
                apply(rsample, 2, function(col) sum(col < -scopes$q) >= sum(scopes$hatLC) ))
  IVkabs <- mean(apply(rsample, 2, function(col) sum(abs(col) > scopes$q) >= sum(scopes$hatUC) + sum(scopes$hatLC) ))
  IV0kN <- mean(apply(rsample, 2, function(col) sum(abs(col) > scopes$kN) >= sum(scopes$hatUC) + sum(scopes$hatLC)  ))
  IV  <- vapply(1:10, function(k) mean(apply(rsample, 2, function(col) sum(abs(col) > scopes$q) >= k )), FUN.VALUE = 0.1)
  IVo <- t(vapply( Tnull, function(t)
        vapply(1:10, function(k) mean(apply(rsample, 2, function(col) sum(abs(col) > t) >= k )), FUN.VALUE = 0.1),
        FUN.VALUE = rep(0.1, 10)))
  colnames(IVo) <- 1:10

  IVo2 <- t(vapply( Tnull, function(t){
    lx = length(scopes$x)
#    mm = ifelse(hatm0 != 0, hatm0, 1)
#    vapply(1:10, function(k) mean(apply(rsample[1:mm, ], 2,
#                                        function(col) sum(abs(col) > t) >= k )), FUN.VALUE = 0.1)},
    vapply(1:10, function(k) mean(apply(rsample[1:(lx - detectN + 1), ], 2,
                      function(col) sum(abs(col) > t) >= k )), FUN.VALUE = 0.1)},
    FUN.VALUE = rep(0.1, 10)))
  IVo3 <- vapply( (0:(length(scopes$x)-2)), function(k)
                  mean(colSums(abs(rsample[1:(length(scopes$x) - k), ] > scopes$q)) >= 1),
                  FUN.VALUE = 0.1)
  names(IVo3) <- length(scopes$x) - 0:(length(scopes$x)-2)
  # return the results
  if(!is.null(scopes$kN)){
    # Global probability of rejection of Gamma(mu) in Gamma(C)
    # if hatmu1C is non-empty
    IV0 = 1 - cdf_abs(scopes$kN)^length(scopes$x)

    return(list(IV0 = IV0, IV0kN = IV0kN, IV = IV,  IVo = IVo, IVo2 = IVo2, IVo3 = IVo3, IVk = IVk, IVkabs = IVkabs, IVobs = IVobs,
                IVloc = IVloc, hatm0 = hatm0))
  }else{
    return(list(IV = IV, IVo = IVo, IVo2 = IVo2, IVo3 = IVo3, IVk = IVk, IVkabs = IVkabs, IVobs = IVobs, IVloc = IVloc, hatm0 = hatm0))
  }
}
