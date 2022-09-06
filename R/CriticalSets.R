#------------------------------------------------------------------------------#
#                                                                              #
#     Functions to estimate critical sets for CoPE sets or Dette's tests       #                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#
#------------------------------------------------------------------------------#
# Developer notes:
#
#------------------------------------------------------------------------------#
#' Estimates the .
#'
#' @param hatmu function to be thresholded to estimated the set SC. Usually it is
#' an estimator of the signal
#' @param C T x 2 matrix containing in the first column the lower bound for the
#' tube and in the second column the upper bound.
#' @param hatsigma standard deviation used to estimate the set SC. Usually it is
#' the pointwise standard deviation
#' @param tN is inverse of the rate of the fCLT, i.e., usually 1/sqrt(N)
#' @param  kN is the thickening of the set (Default log(tN^-2) / 5)
#' @param  type if CoPE sets are computed either "selection" or
#' "extraction"(Default ). Otherwise if it is estimated for a statistical test
#' using CoPE sets it is "rel" (relevant tube test), "out" (outside relevant tube
#'  test), "eqv" (bio-equivalence test)
#' @return list with elements
#'  \itemize{
#'   \item hatmean pointwise sample mean
#'   \item scb list containing the upper and lower bounds of the simultaneous confidence band
#'   \item level targeted covering probability
#'   \item q quantile of the maximum of the residual field
#' }
#' @export
PreimageC <- function(hatmu,
                      C,
                      hatsigma,
                      tN,
                      kN = log(tN^-2) / 5,
                      type = "extraction"){
  # Get the dimension of the
  dC = dim(C)
  if(type %in% c("rel", "out", "eqv")){
    # Construct the threshold of the CoPE set
    if(type == "rel"){
      Delta = min(min(abs(C[,1] - hatmu)), min(abs(hatmu - C[,dC[2]])))
      cm = C[,1] + Delta
      cp = C[,dC[2]] - Delta

      updn <- c(diff(sign(hatmu-C[,1])), 0)
      ix <- which(updn != 0)
      ixm = vapply(ix, function(v){ ifelse( which.min(c(abs(hatmu[v]-C[v,1]),
                                                        abs(hatmu[v+1]-C[v+1,1]))) == 1,
                                            v, v+1)}, FUN.VALUE = 1 )
      updn <- c(diff(sign(hatmu-C[,dC[2]])), 0)
      ix <- which(updn != 0)
      ixp = vapply(ix, function(v){ ifelse( which.min(c(abs(hatmu[v]-C[v,dC[2]]),
                                                        abs(hatmu[v+1]-C[v+1,dC[2]]))) == 1,
                                            v, v+1)}, FUN.VALUE = 1 )

    }else if(type == "out"){
      Delta = min(min(abs(C[,1] - hatmu)), min(abs(hatmu - C[,dC[2]])))
      cm = C[,1] - Delta
      cp = C[,dC[2]] + Delta
      ixm <- ixp <- NULL
    }else if(type == "eqv"){
      Delta = max(max(C[,1] - hatmu), max(hatmu - C[,dC[2]]))
      cm = C[,1] - Delta
      cp = C[,dC[2]] + Delta
      ixm <- ixp <- NULL
    }

    # Get the estimates of the set
    S_cm = ( cm - tN * kN * hatsigma <= hatmu) & (hatmu <= cm + tN * kN * hatsigma)
    S_cp = ( cp - tN * kN * hatsigma <= hatmu) & (hatmu <= cp + tN * kN * hatsigma)

    S_cm[ixm] <- TRUE
    S_cp[ixp] <- TRUE
    # Output the set S_p and S_m
    list(
      S_m = S_cm,
      S_p = S_cp,
      c_m = cm,
      c_p = cp,
      bd  = tN * kN * hatsigma
    )
  }else if(type == "extraction"){
    M = (C - tN * kN * hatsigma <= hatmu) & (hatmu <= C + tN * kN * hatsigma)
    out = list(
      S_p = M[,1],
      S_m = M[,dC[2]]
    )
    return( out )
    # return( apply(M, 1, any) )
  }else if(type == "selection"){
    return( (C[,1] - tN * kN * hatsigma <= hatmu) & (hatmu <= C[,dC[2]] + tN * kN * hatsigma) )
  }
}

#' Estimates the critical set for the Dette papers.
#' @inheritParams PreimageC
#' @return list with elements
#'  \itemize{
#'   \item hatmean pointwise sample mean
#'   \item scb list containing the upper and lower bounds of the simultaneous confidence band
#'   \item level targeted covering probability
#'   \item q quantile of the maximum of the residual field
#' }
#' @export
# hat Extremal sets E_pm
hatE_pm <- function(hatmu, C, hatsigma, tN, kN = log(tN^-2) / 5){
  T_hatmu = max(max(-hatmu + C[,1]), max(hatmu - C[,2]))

  out = list(
    S_p =  hatmu - C[,2] >= T_hatmu - tN * kN * hatsigma,
    S_m = -hatmu + C[,1] >= T_hatmu - tN * kN * hatsigma
  )
  return(out)
}
