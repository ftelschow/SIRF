#------------------------------------------------------------------------------#
#                                                                              #
#     Functions to compute CoPE sets for functional data
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#
#
#------------------------------------------------------------------------------#
# Developer notes:
#
#------------------------------------------------------------------------------#
#' This functions computes an lower confidence excursion set
#'
#' @inheritParams PreimageC
#' @param c vector with the same length as hatmu containing the thresholding
#' function
#' @param q quantile for the CoPE set
#' @param inclusion if "inequal" the lower excursion set is open, if "equal"
#' the lowerer excursion set is closed.
#' @return  if is.null(mu) a boolean vector indicating which positions are
#' in the lower confidence excursion set. Otherwise a list with elements
#'  \itemize{
#'   \item hatUc boolean vector indicating which positions are
#' in the lower confidence excursion set
#'   \item contain boolean indicating whether the lower confidence set is contained in
#'   the true lower excursion set
#' }
#' @export
hatLc <- function(c, hatmu, hatsigma, tN, q, mu = NULL, inclusion = "inequal" ){
  if(inclusion == "inequal"){
    hatL = hatmu < c - q * tN * hatsigma
    if(!is.null(mu)){
      Lc = mu < c
      contain  = (hatL | Lc) == Lc
      return(list(hatLc = hatL, contain = contain))
    }else{
      return(hatL)
    }
  }else{
    hatU = hatUc(c = c, hatmu = hatmu, hatsigma = hatsigma,
                 tN = tN, q = q, mu = mu, inclusion = "inequal" )
    if(!is.null(mu)){
      return(list(hatLc = !hatU$hatUc, contain = hatU$contain))
    }else{
      return(!hatU)
    }
  }
}

#' This functions computes an upper confidence excursion set
#'
#' @inheritParams PreimageC
#' @param c vector with the same length as hatmu containing the thresholding function
#' @param q quantile for the CoPE set
#' @param inclusion if "inequal" the upper excursion set is open, if "equal" the upper
#' excursion set is closed.
#' @return  if is.null(mu) a boolean vector indicating which positions are
#' in the upper confidence excursion set. Otherwise a list with elements
#'  \itemize{
#'   \item hatUc boolean vector indicating which positions are
#' in the upper confidence excursion set
#'   \item contain boolean indicating whether the upper confidence set is contained in
#'   the true upper excursion set
#' }
#' @export
hatUc <- function(c, hatmu, hatsigma, tN, q, mu = NULL, inclusion = "inequal" ){
  if(inclusion == "inequal"){
    hatU = hatmu > c + q * tN * hatsigma
    if(!is.null(mu)){
      Uc = mu > c
      contain  = (hatU | Uc) == Uc
      return(list(hatUc = hatU, contain = contain))
    }else{
      return(hatU)
    }
  }else{
    hatL = hatLc(c = c, hatmu = hatmu, hatsigma = hatsigma,
                 tN = tN, q = q, mu = mu, inclusion = "inequal" )
    if(!is.null(mu)){
      return(list(hatUc = !hatL$hatLc, contain = hatL$contain))
    }else{
      return(!hatL)
    }
  }
}

#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams PreimageC
#' @param alpha numeric between 0 and 1. It will produce (1-alpha) SCoPES.
#' @param inclusionL determines whether the lower excursion set uses < c or <=c.
#' @param inclusionU determines whether the upper excursion set uses > c or >=c.
#' @return Standard error under the assumption the data is Gaussian
#' @export
SCoPES <- function(alpha, C, hatmu, hatsigma, tN, method = "extraction", R = NULL,
                   inclusionL = "inequal", inclusionU = "inequal",
                   mu = NULL){
  if(method %in% c("SCB", "selection", "extraction", "lrelevance", "relevance",
                   "lequivalence", "equivalence") && !is.null(R) ){
    method <- list(name      = "mboot",
                   Boottype  = "t",
                   weights   = "rademacher",
                   Mboots    = 2e4,
                   R         = R,
                   SCoPEtype = method,
                   mu1Cest   = "thickening",
                   kN        = log(tN^-2) / 5
                   )
  }

  # Get the number of levels and make C a column matrix if
  # it is just one function.
  if(is.null(dim(C))){
    num_levels = 1
    C = t(t(C))
  }else{
    num_levels = dim(C)[2]
  }

  # Estimate the preimage using the true set or the thickened set
  if(is.null(method$mu1Cest)){
    hatmu1C =  PreimageC(C = C, hatmu = mu, hatsigma = rep(0, length(hatmu)),
                         tN = tN, kN = 0, method = method$SCoPEtype)
  }else if(method$mu1Cest == "thickening"){
    hatmu1C =  PreimageC(C = C, hatmu = hatmu, hatsigma = hatsigma,
                         tN = tN, kN = method$kN, method = method$SCoPEtype)
  }

  # Estimate the quantile for the SCoPES
  if(method$name == "mboot"){
    q = MultiplierBootstrapSplit(alpha   = alpha,
                                 R       = method$R,
                                 minus   = hatmu1C$minus,
                                 plus    = hatmu1C$plus,
                                 Mboots  = method$Mboots,
                                 method  = method$Boottype,
                                 weights = method$weights)$q
  }else if(method$name == "gauss"){
    q = maxGauss_quantile(1 - alpha, hatmu1C)
  }else if(method$name == "t"){
    q = maxT_quantile(1 - alpha, hatmu1C, df = method$df)
  }

  # Initialize matrizes for the CoPE sets and for the rejections
  hatLC    <- matrix(FALSE, length(hatmu), num_levels)
  hatUC    <- matrix(FALSE, length(hatmu), num_levels)
  containL <- matrix(FALSE, length(hatmu), num_levels)
  containU <- matrix(FALSE, length(hatmu), num_levels)

  for(k in 1:num_levels){
    hatLc <- hatLc(C[, k], hatmu, hatsigma, tN, q, mu = mu, inclusion = inclusionL )
    hatUc <- hatUc(C[, k], hatmu, hatsigma, tN, q, mu = mu, inclusion = inclusionU )

    if(!is.null(mu)){
      hatLC[, k]    <- hatLc$hatLc
      containL[, k] <- hatLc$contain
      hatUC[, k]    <- hatUc$hatUc
      containU[, k] <- hatUc$contain

    }else{
      hatLC[, k] <- hatLc
      hatUC[, k] <- hatUc
    }
  }

  if(!is.null(mu)){
    return(list(hatLC = hatLC, hatUC = hatUC, q = q, hatmu1C = hatmu1C,
                containL = containL, containU = containU))
  }else{
    return(list(hatLC = hatLC, hatUC = hatUC, q = q, hatmu1C = hatmu1C))
  }
}
