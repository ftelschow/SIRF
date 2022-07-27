#------------------------------------------------------------------------------#
#                                                                              #
#     Functions to compute CoPE sets for functional data
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#  - scb_moments
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline included
# - add method description
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
PreimageC <- function(hatmu, C, hatsigma, tN, kN = log(tN^-2) / 5, type = "extraction"){
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


#' This functions computes inner CoPE set
#'
#' @inheritParams PreimageC
#' @param c vector with the same length as hatmu contianing the thresholding function
#' @param q quantile for the CoPE set
#' @param inclusion if "inequal" the excursion set is open, if "equal" the excursion set
#' is closed.
#' @return Standard error under the assumption the data is Gaussian
#' @export
inner_CoPEset <- function(hatmu, c, q, tN, hatsigma, inclusion = "inequal", mu = NULL ){
  if(inclusion == "equal"){
    innerSet = hatmu >= c + q * tN * hatsigma
    if(!is.null(mu)){
      Ac = mu >= c
      innerContain  = all((innerSet | Ac) == Ac)
    }else{
      innerContain = NaN
    }
  }else{
    innerSet = hatmu > c + q * tN * hatsigma
    if(!is.null(mu)){
      Ac = mu > c
      innerContain  = all((innerSet | Ac) == Ac)
    }else{
      innerContain = NaN
    }
  }
  list(innerSet = innerSet, contain = innerContain)
}

#' This functions computes inner CoPE set
#'
#' @inheritParams PreimageC
#' @param c vector with the same length as hatmu contianing the thresholding function
#' @param q quantile for the CoPE set
#' @param inclusion if "inequal" the excursion set is open, if "equal" the excursion set
#' is closed.
#' @return Standard error under the assumption the data is Gaussian
#' @export
outer_CoPEset <- function(hatmu, c, q, tN, hatsigma, inclusion = "inequal", mu = NULL ){
  if(inclusion == "equal"){
    outerSet = hatmu >= c - q * tN * hatsigma
    if(!is.null(mu)){
      Ac = mu >= c
      outerContain  = all((outerSet | Ac) == outerSet)
    }else{
      outerContain = NaN
    }
  }else{
    outerSet = hatmu > c - q * tN * hatsigma
    if(!is.null(mu)){
      Ac = mu > c
      outerContain  = all((outerSet | Ac) == outerSet)
    }else{
      outerContain = NaN
    }
  }
  list(outerSet = outerSet, contain = outerContain)
}

#' This functions computes the CoPE sets
#'
#' @inheritParams PreimageC
#' @param inclusion determines whether the CoPE sets are closed or open.
#' The usually correct choice is "inequalequal" (default), i.e., the inner
#' CoPE set is open, the outer is closed. Other options are "equalinequal",
#' "equal", "inequal".
#' @return Standard error under the assumption the data is Gaussian
#' @export
CoPEsets <- function(hatmu, C, q, tN, hatsigma, inclusion = "inequalequal", mu = NULL, eps_correct = T, subI = NULL){
  # Get the number of levels
  if(is.null(dim(C))){
    num_levels = 1
    C = t(t(C))
  }else{
    num_levels = dim(C)[2]
  }

  if(is.null(subI)){
    subI = list()
    for(k in 2:length(fair)){
      subI[[k-1]] = which( x  >= fair[k-1] & x  <= fair[k] )
    }
  }



  # Initialize matrizes for the CoPE sets and for the rejections
  innerSet <- outerSet <- matrix(NaN, length(hatmu), num_levels)
  contain  <- matrix(NaN, 2, num_levels)
  row.names(contain) = c("inner", "outer")

  for(k in 1:num_levels){
    innerContain = TRUE
    outerContain = TRUE
    if(eps_correct){
      if(k == 1)
        epsVec = c(1e-8, 0)
      else if(k==num_levels){
        epsVec = c(-1e-8, 0)
      }else{
        epsVec = c(-1e-8, 1e-8, 0)
      }
    }else{
      epsVec = 0
    }

    for(eps in epsVec){
      c = C[, k] + eps
      if(inclusion == "equal"){
        innerSet.c = inner_CoPEset(hatmu, c, q, tN, hatsigma, inclusion = "equal", mu = mu )
        outerSet.c = outer_CoPEset(hatmu, c, q, tN, hatsigma, inclusion = "equal", mu = mu )
      }else if(inclusion == "inequal"){
        innerSet.c = inner_CoPEset(hatmu, c, q, tN, hatsigma, inclusion = "inequal", mu = mu )
        outerSet.c = outer_CoPEset(hatmu, c, q, tN, hatsigma, inclusion = "inequal", mu = mu )
      }else if(inclusion == "equalinequal"){
        innerSet.c = inner_CoPEset(hatmu, c, q, tN, hatsigma, inclusion = "equal", mu = mu )
        outerSet.c = outer_CoPEset(hatmu, c, q, tN, hatsigma, inclusion = "inequal", mu = mu )
      }else if(inclusion == "inequalequal"){
        innerSet.c = inner_CoPEset(hatmu, c, q, tN, hatsigma, inclusion = "inequal", mu = mu )
        outerSet.c = outer_CoPEset(hatmu, c, q, tN, hatsigma, inclusion = "equal", mu = mu )
      }
      innerContain = innerContain & innerSet.c$contain
      outerContain = outerContain & outerSet.c$contain
    }

    innerSet[, k] = innerSet.c$innerSet
    outerSet[, k] = outerSet.c$outerSet

    contain[1, k] = innerContain
    contain[2, k] = outerContain
  }

  list(innerSet = innerSet, outerSet = outerSet, contain = contain)
}

#' This functions enables different inference for 1D fields using
#' CoPE sets
#'
#' @inheritParams PreimageC
#' @param x The considered transformation of the moments
#' @param scb Vector grid for evaluation of the smoothed field
#' @param scb2 Vector grid for evaluation of the smoothed field
#' @param true Vector grid for evaluation of the smoothed field
#' @param title Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
CoPE_inference <- function(hatmu, hatsigma, R,
                               type,
                               C,
                               alpha  = 0.05,
                               tN,
                               kN,
                               Mboots = 1e4,
                               mu = NULL,
                               fair = NULL,
                               fair.type = "linear",
                               niter  = 10,
                               subI = NULL,
                               print.coverage = TRUE
){
  #
  if(kN == 0 & !is.null(mu)){
    estmu = mu
    method = "regular"
    weights = "gaussian"
  }else{
    estmu = hatmu
    method = "t"
    weights = "rademacher"
  }

  if(is.null(subI) && !is.null(fair)){
    subI = list()
    for(k in 2:length(fair)){
      subI[[k-1]] = which( x  >= fair[k-1] & x  <= fair[k] )
    }
  }else{
    subI = NULL
  }


  # Estimate the quantile. This requires estimation of the critical set
  if(type == "extraction"){
    SC = PreimageC(hatmu = estmu, C = C, hatsigma = hatsigma,
                   tN = tN, kN = kN, type = "rel")
    Splus  = SC$S_p
    Sminus = SC$S_m

    # Get the unfair threshold
    q = MultiplierBootstrapSplit(R = R, Splus = Splus,
                                 Sminus = Sminus,
                                 alpha  = alpha,
                                 Mboots = Mboots,
                                 method  = method,
                                 weights = weights)

    # Get the optimized fair piecewise linear threshold
    if(!is.null(fair)){
      test = OptimizeFairThreshold1D(q$samples, x, fair,
                                     fair.type = fair.type,
                                     Splus  = Splus,
                                     Sminus = Sminus,
                                     alpha  = alpha,
                                     niter  = niter,
                                     print.coverage = print.coverage,
                                     subI = subI)
      q$q = test$q
      q$EmpRejections = test$EmpRejections
    }

    # ColorVector of estimated set
    colVecSet = rep("black", length(x))
    colVecSet[SC$S_p] = "red"
    colVecSet[SC$S_m] = "blue"

  }else if(type == "selection"){
    SC = PreimageC(hatmu = estmu, C = C, hatsigma = hatsigma,
                   tN = tN, kN = kN, type = "selection")

    if(is.null(fair)){
    q = MultiplierBootstrap(R = R[SC,],
                            alpha  = alpha,
                            Mboots = Mboots,
                            method  = method,
                            weights = weights)
   }else{
      q = MultiplierBootstrap(R = R,
                             alpha  = alpha,
                             Mboots = Mboots,
                             method  = method,
                             weights = weights)
      test = OptimizeFairThreshold1D(q$samples, x, fair,
                                     fair.type = fair.type,
                                     Splus  = SC,
                                     Sminus = SC,
                                     alpha  = alpha,
                                     niter  = niter,
                                     print.coverage = print.coverage,
                                     subI = subI )
      q$q = test$q
      q$EmpRejections = test$EmpRejections
    }
    # ColorVector of estimated set
    colVecSet = rep("black", length(x))
    colVecSet[SC] = "red"

  }else if(type == "classical"){
    SC = PreimageC(hatmu = estmu, C = C, hatsigma = hatsigma,
                   tN = tN, kN = kN, type = "rel")
    SCC = SC$S_p | SC$S_m
    if(sum(SCC) == 1){
      # Get the threshold
      q = MultiplierBootstrap(R = t(R[SCC,]),
                              Mboots = Mboots,
                              alpha  = alpha,
                              method  = method,
                              weights = weights)

    }else if(sum(SCC) > 1){
      if(is.null(fair)){# unfair threshold
        q = MultiplierBootstrap(R = t(t(R[SCC,])),
                                Mboots = Mboots,
                                alpha  = alpha,
                                method  = method,
                                weights = weights)

      }else{# fair threshold
        q = MultiplierBootstrap(R = R,
                                alpha  = alpha,
                                Mboots = Mboots,
                                method  = method,
                                weights = weights)
        test = OptimizeFairThreshold1D(q$samples, x, fair,
                                       fair.type = fair.type,
                                       Splus  = SCC,
                                       Sminus = SCC,
                                       alpha  = alpha,
                                       niter  = niter,
                                       print.coverage = print.coverage,
                                       subI = subI)
        q$q = test$q
        q$EmpRejections = test$EmpRejections
      }
    }else{
      q = list()
      q$q = 0
    }
    # ColorVector of estimated set
    colVecSet = rep("black", length(x))
    colVecSet[SC$S_m | SC$S_p] = "red"
  }else if(type == "SCB"){
    # Get the unfair threshold
    q = MultiplierBootstrap(R = R,
                            Mboots = Mboots,
                            alpha = alpha,
                            method  = method,
                            weights = weights)

    # Get the optimized fair piecewise linear threshold
    if(!is.null(fair)){
      test = OptimizeFairThreshold1D(q$samples, x, fair,
                                     fair.type = fair.type,
                                     alpha  = alpha,
                                     niter  = niter,
                                     print.coverage = print.coverage,
                                     subI = subI)
      q$q = test$q
      q$EmpRejections = test$EmpRejections
    }

    # ColorVector of estimated set
    colVecSet = rep("red", length(x))
  }

  # Modify C such that it not automatically rejects, if q is -Inf because there
  # are no elements of S_C
  if(!is.null(fair)){
    C[which(is.infinite(q$q)), ] = -Inf
  }

  # Get the CoPE sets
  coPE <- CoPEsets(hatmu = mY,
                   C = C,
                   q = q$q,
                   hatsigma = sdY,
                   tN = tN,
                   mu = mu)

  colVecCoPE = rep("black", length(x))
  colVecCoPE[!!coPE$innerSet[,dim(coPE$innerSet)[2]]] = "red"
  colVecCoPE[!coPE$outerSet[,1]]  = "blue"
  Reject = list()
  Reject$low = which(!coPE$outerSet[,1])
  Reject$up  = which(!!coPE$innerSet[,dim(coPE$innerSet)[2]])

  if(is.null(fair)){
    list(Reject  = Reject,
         CoPEcol = colVecCoPE,
         SetCol  = colVecSet,
         CoPE    = coPE,
         q = q,
         CoPE_bands = hatmu + tN * q$q * cbind(-hatsigma, hatsigma))
  }else{
    list(Reject  = Reject,
         CoPEcol = colVecCoPE,
         SetCol  = colVecSet,
         CoPE    = coPE,
         q       = q,
         C       = C,
         loc_alpha  = test$loc_alpha,
         CoPE_bands = hatmu + tN *  cbind(-q$q * hatsigma, q$q * hatsigma) )
  }

}

#' This functions enables different inference for 1D fields using
#' CoPE sets
#'
#' @inheritParams PreimageC
#' @return Standard error under the assumption the data is Gaussian
#' @export
Dette_inference <- function(hatmu, hatsigma, R,
                           C,
                           alpha  = 0.05,
                           tN,
                           kN,
                           Mboots = 1e4,
                           mu = NULL
){
  # Get the critical set of Dette
  SC = hatE_pm(hatmu = hatmu, C = C, hatsigma = hatsigma, tN = tN, kN = kN)
  Splus  = SC$S_p
  Sminus = SC$S_m

  # Bootstrap the quantile
  q = MultiplierBootstrapSplit(R = R, Splus = Splus,
                               Sminus = Sminus,
                               alpha  = alpha,
                               Mboots = Mboots)

   # ColorVector of estimated set
   colVecSet = rep("black", length(x))
   colVecSet[SC$S_p] = "red"
   colVecSet[SC$S_m] = "blue"

  # Dette's Test statistic
  T_hatmu = max(max(-hatmu + C[,1]), max(hatmu - C[,2]))

  # Compute the CoPE sets for completeness
  coPE <- CoPEsets(hatmu = mY,
                   C = C,
                   q = q$q,
                   hatsigma = sdY,
                   tN = tN,
                   mu = mu)

  colVecCoPE = rep("black", length(x))
  colVecCoPE[!!coPE$innerSet[,2]] = "red"
  colVecCoPE[!coPE$outerSet[,1]]  = "blue"

  # Give the output of the bio-equivalence and the relevance test
  Reject = c("accept", "accept")
  names(Reject) <- c("equivalence", "relevance")
  if(all(colVecCoPE == "black")){
    Reject[1] <- "reject"
  }else{
    Reject[2] <- "reject"
  }

  list(Reject  = Reject,
       CoPEcol = colVecCoPE,
       SetCol  = colVecSet,
       CoPE    = coPE,
       q = q)
}
