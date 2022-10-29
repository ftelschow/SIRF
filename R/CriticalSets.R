#------------------------------------------------------------------------------#
#                                                                              #
#     Functions to estimate critical sets for CoPE sets or Dette's tests       #
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#
#------------------------------------------------------------------------------#
# Developer notes:
#
#------------------------------------------------------------------------------#
#' Estimates the preimage of the set indicated by the matrix C using
#' the thickening estimator proposed in Telschow et al (2022)
#' "SCoPES: A versatile framework of simultaneous inference".
#'
#' @param C T x nlvl matrix containing nlvl functions for which the preimage
#'          should be computed. At the moment only the first and last function
#'          matter, i.e., they define a tube.
#' @param hatmu function to be thresholded to estimated the set mu1_C.
#' Usually it is an estimator of the target function mu
#' @param hatsigma scaling used to estimate the set mu1_C. Usually it is
#' the pointwise standard deviation og hatmu.
#' @param tN is inverse of the rate of the fCLT, i.e., usually 1/sqrt(N)
#' @param kN is the thickening of the set (Default log(tN^-2) / 5)
#' @param method if CoPE sets are computed either "selection" or
#' "extraction"(Default ).
#'  Otherwise if it is estimated for a statistical test
#'  using CoPE sets the options are "lrelevance", "relevance",
#'  "lequivalence", "equivalence"
#' @return list with elements
#'  \itemize{
#'   \item minus boolean vector indicating which locations are in the set
#'         mu1_C^-
#'   \item plus boolean vector indicating which locations are in the set
#'         mu1_C^+
#' }
#' @export
PreimageC <- function(C,
                      hatmu,
                      hatsigma,
                      tN,
                      kN     = log(tN^-2) / 5,
                      method = "extraction"){
  # Get the dimension of the preimage matrix
  dC = dim(C)

  # Delta for moving the boundary for testing statistical hypotheses
  Delta = 0
  # Correct Delta for testing
  if(method %in%  c("lrelevance", "lequivalence")){
    Delta = min(min(abs(C[, 1] - hatmu)),
                min(abs(hatmu - C[, dC[2]])))
    if(method == "lrelevance"){
      C[, 1]     = C[, 1]     + Delta
      C[, dC[2]] = C[, dC[2]] - Delta
    }else{
      C[, 1]           = C[, 1]     - Delta
      C[, dC[2]]       = C[, dC[2]] + Delta
      C[, c(1, dC[2])] = C[, c(dC[2], 1)]
    }
  }else if(method %in%  c("equivalence", "relevance")){
    Delta = max(max(C[, 1] - hatmu),
                max(hatmu - C[, dC[2]]))
    C[, 1]     = C[, 1]     - Delta
    C[, dC[2]] = C[, dC[2]] + Delta
  }

  # Get the estimate of the preimages from the estimated
  # preimages of C for the different methods
  if(method %in% c("extraction", "lrelevance", "relevance",
                   "lequivalence", "equivalence")){
    # Get the estimate of the preimage for each col \in C
    M = (C[, c(1, dC[2])] - hatmu <= tN * kN * hatsigma) &
              (hatmu - C[, c(1, dC[2])] <= tN * kN * hatsigma)

    hatmu1_C = list(minus = M[, 1], plus  = M[, 2])

  }else if(method == "selection"){
    # Get the estimate of the preimage for each col \in C
    M = (C[, 1] - hatmu <= tN * kN * hatsigma) &
              (hatmu - C[, dC[2]] <= tN * kN * hatsigma)
    hatmu1_C = list(minus = M, plus  = M)
  }else{
    stop("The input method must be from the following list of
         options: extraction, selection, lrelevance, relevance,
         lequivalence, equivalence.")
  }
  # Return the estimate of the preimage of C
  return(hatmu1_C)
}
