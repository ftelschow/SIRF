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
      detect  <- NA * (1:length(hatmu))
      detect[hatL] <- contain[hatL]
      return(list(hatLc = hatL, contain = contain, detect = detect))
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
      contain <- (hatU | Uc) == Uc
      detect  <- NA * (1:length(hatmu))
      detect[hatU] <- contain[hatU]
      return(list(hatUc = hatU, contain = contain, detect = detect))
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
#' @param method this can be either a string specifying the method, i.e.,
#'  "selection", "extraction"(Default ), lrelevance", "relevance",
#'  "lequivalence", "equivalence". The latter for are the hypothesis tests
#'  described in Telschow et al 2022. In this case the quantile q is obtained
#'  by a bootstrap on the thickened set.
#'  It can also be a method specifying the estimation of the quantile $q$.
#' @param R array containing as columns residuals which can be used to bootstrap
#' the quantile q. This input is mandatory if "method" is only a string.
#' @param inclusion list with elements
#'  \itemize{
#'   \item L string either "inequal" or "equal" determining whether the lower
#'   excursion set uses < c or <=c.
#'   \item U string either "inequal" or "equal" determining whether the upper
#'   excursion set uses > c or >=c.
#' }
#' @return Standard error under the assumption the data is Gaussian
#' @export
SCoPES <- function(alpha, C, x, hatmu, hatsigma, tN, method = "extraction", R = NULL,
                   I = NULL, inclusion = list(L = "inequal", U = "inequal"),
                   mu = NULL){
  #-----------------------------------------------------------------------------
  # Get standard values and catch wrong inputs

  # Get the number of levels and make C a column matrix if
  # it is just one function.
  if(is.null(dim(C))){
    num_levels = 1
    C = t(t(C))
  }else{
    num_levels = dim(C)[2]
  }

  # Dimension of the tube defining matrix
  dC <- dim(C)
  rownames(C) <- round(x,3)
  colnames(C) <- 1:dC[2]
  Cnames <- list(x = round(x,3), c = 1:dC[2])

  # Get standard values for method, if not specified
  if(method %in% c("SCB", "selection", "extraction", "lrelevance", "relevance",
                   "lequivalence", "equivalence") && !is.null(R) ){
    method <- list(name      = "mboot",
                   Boottype  = "t",
                   weights   = "rademacher",
                   Mboots    = 2e4,
                   R         = R,
                   SCoPEStype = method,
                   mu1Cest   = "thickening",
                   kN        = log(tN^-2) / 5
                   )
  }
  # Get standard values for the indices of C^\pm if not specified
  if(is.null(I)){
    if(method$SCoPEStype %in% c( "extraction", "relevance", "lrelevance",
                         "equivalence", "lequivalence")){
      Iminus <- seq(1, dC[2], 2)
      Iplus  <- seq(2, dC[2], 2)
    }else{
      Iminus <- 1:dC[2]
      Iplus  <- 1:dC[2]
    }
    I = list(minus = Iminus, plus = Iplus)
  }else if(is.list(I)){
    if(length(I) == 2){
      Iminus = I$minus
      Iminus = I$plus
    }else{
      stop("I must be a list with two entries indicating which columns of
           C belong to C^- and which to C^+")
    }
  }else{
    stop("I must be a list with two entries indicating which columns of
           C belong to C^- and which to C^+")
  }

  # Ensure that for the tests C defines a single tube
  if(method$SCoPEStype %in% c("relevance", "lrelevance",
                              "equivalence", "lequivalence")){
    if(dC[2] > 2){
      stop("For 'relevance'', 'lrelevance', 'equivalence', 'lequivalence' C must
           have two columns defining the considered band.")
    }
  }else if(method$SCoPEStype == "classical"){
    if(num_levels > 1){
      stop("For 'classical' SCoPES C must be a vector or column matrix with 1
           column.")
    }
  }

  #-----------------------------------------------------------------------------
  # Main code of the algorithm

  # Estimate the preimage using the true set or the thickened set
  if(method$SCoPEStype != "SCB"){
    if(is.null(method$mu1Cest)){
      hatmu1C =  PreimageC(C = C, hatmu = mu, hatsigma = rep(0, length(hatmu)),
                           tN = tN, kN = 0, method = method$SCoPEStype)
    }else if(is.list(method$mu1Cest)){
      hatmu1C = method$mu1Cest
    }else if(method$mu1Cest == "thickening"){
      hatmu1C =  PreimageC(C = C, hatmu = hatmu, hatsigma = hatsigma,
                           tN = tN, kN = method$kN, method = method$SCoPEStype)
    }else if(method$mu1Cest == "m0"){
      hatmu1C = method$m0
    }
  }else{
    hatmu1C = list(minus = rep(TRUE, length(hatmu)), plus = rep(TRUE, length(hatmu)))
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
    if(is.infinite(q)){
      q = 0
    }
  }else if(method$name == "gauss"){
    q = maxGauss_quantile(p = 1 - alpha, muC = hatmu1C)
  }else if(method$name == "t"){
    q = maxT_quantile(p = 1 - alpha, muC = hatmu1C, df = method$df)
  }

  # Initialize matrizes for the CoPE sets and for the rejections
  hatLC        <- t(t(matrix(FALSE, length(hatmu), length(Iminus))))
  hatUC        <- t(t(matrix(FALSE, length(hatmu), length(Iplus))))
  Lcontain_loc <- t(t(matrix(FALSE, length(hatmu), length(Iminus))))
  Ucontain_loc <- t(t(matrix(FALSE, length(hatmu), length(Iplus))))
  DetectL  <- t(t(matrix(NA, length(hatmu), length(Iminus))))
  DetectU  <- t(t(matrix(NA, length(hatmu), length(Iplus))))

  k_minus <- k_plus <- 0
  for(k in 1:num_levels){
    # Compute the lower and upper excursion confidence set, if required
    if(k %in% Iminus){
      hatL <- hatLc(C[, k], hatmu, hatsigma, tN, q, mu = mu, inclusion = inclusion$L )
      k_minus <- k_minus + 1
      # Fill the hatLC and estUc matrix and save where correct inclusions appear
      if(!is.null(mu)){
        hatLC[, k_minus]        <- hatL$hatLc
        Lcontain_loc[, k_minus] <- hatL$contain
        DetectL[, k_minus]      <- hatL$detect
      }else{
        hatLC[, k_minus] <- hatL
      }
    }
    if(k %in% Iplus){
      hatU <- hatUc(C[, k], hatmu, hatsigma, tN, q, mu = mu, inclusion = inclusion$U )
      k_plus <- k_plus + 1
      # Fill the hatLC and estUc matrix and save where correct inclusions appear
      if(!is.null(mu)){
        hatUC[, k_plus]         <- hatU$hatUc
        Ucontain_loc[, k_plus]  <- hatU$contain
        DetectU[, k_plus]       <- hatU$detect
      }else{
        hatUC[, k_plus]  <- hatU
      }
    }
  }

  if(!is.null(mu)){
    # Compute the coverage properties
    if(method$SCoPEStype %in% c("relevance", "lrelevance",
                                  "equivalence", "lequivalence")){
      Lcontain_loc <- Lcontain_loc
      Ucontain_loc <- Ucontain_loc
    }else if(method$SCoPEStype %in% c("extraction", "classical")){
      Lcontain_loc <- apply(t(t(Lcontain_loc)), 1, all)
      Ucontain_loc <- apply(t(t(Ucontain_loc)), 1, all)
    }else if(method$SCoPEStype %in% c("selection", "SCB")){
      Lcontain_loc <- apply(t(t(Lcontain_loc)), 1, all)
      Ucontain_loc <- apply(t(t(Ucontain_loc)), 1, all)
    }else{
      stop("q_method$SCoPEStype must be either 'classical', 'extraction',
            'selection', 'lrelevance', 'relevance', 'lequivalence',
           'equivalence'.")
    }

    # Set the col and row names for the output variables
    rownames(hatLC) <- rownames(hatUC) <- rownames(C)
    colnames(hatLC) <- colnames(C)[Iminus]
    colnames(hatUC) <- colnames(C)[Iplus]

    Lcontain <- all(Lcontain_loc)
    Ucontain <- all(Ucontain_loc)
    contain  <- Lcontain && Ucontain

    t_detect <- sum(t(DetectL), na.rm = T) + sum(t(DetectU), na.rm = T)
    f_detect <- sum(t(DetectL) == 0, na.rm = T) + sum(t(DetectU) == 0, na.rm = T)

    return(list(hatLC = hatLC, hatUC = hatUC, q = q, hatmu1C = hatmu1C, kN = method$kN,
                contain = contain, Lcontain = Lcontain, Ucontain = Ucontain,
                Lcontain_loc = Lcontain_loc, Ucontain_loc = Ucontain_loc,
                LDetect = DetectL, UDetect = DetectU,
                NtrueDetect = t_detect, NfalseDetect = f_detect,
                mu = mu, x = x, tN = tN, hatmu = hatmu, hatsigma = hatsigma,
                C = C, I = I, SCoPEStype = method$SCoPEStype))
  }else{
    return(list(hatLC = hatLC, hatUC = hatUC, q = q, hatmu1C = hatmu1C, kN =  method$kN,
                x = x, tN = tN, hatmu = hatmu, hatsigma = hatsigma,
                C = C, I = I, SCoPEStype = method$SCoPEStype))
  }
}


#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams PreimageC
#' @param alpha numeric between 0 and 1. It will produce (1-alpha) SCoPES.
#' @param method this can be either a string specifying the method, i.e.,
#'  "selection", "extraction"(Default ), lrelevance", "relevance",
#'  "lequivalence", "equivalence". The latter for are the hypothesis tests
#'  described in Telschow et al 2022. In this case the quantile q is obtained
#'  by a bootstrap on the thickened set.
#'  It can also be a method specifying the estimation of the quantile $q$.
#' @param R array containing as columns residuals which can be used to bootstrap
#' the quantile q. This input is mandatory if "method" is only a string.
#' @param inclusion list with elements
#'  \itemize{
#'   \item L string either "inequal" or "equal" determining whether the lower
#'   excursion set uses < c or <=c.
#'   \item U string either "inequal" or "equal" determining whether the upper
#'   excursion set uses > c or >=c.
#' }
#' @return Standard error under the assumption the data is Gaussian
#' @export
plot_SCoPES <- function(scopes, index_C = 1,
                        xlab = '', ylab = '', title = '',
                        mu = NULL, statistic = NULL){
  if(scopes$SCoPEStype == "extraction"){
    index_C <- 2 * index_C - 1
    Delta = 1
  }else{
    Delta = 0
  }
  if(is.null(statistic)){
    y = scopes$hatmu
  }else{
    y = statistic
  }
  # Get the correct
  C = scopes$C[, c(index_C, index_C + Delta)]
  # Get a color vector indicating the "red" (upper excursions) and
  # the blue (lower excursions) set
  colVec <- rep("black", length(scopes$hatmu))
  colVec[scopes$hatLC[, index_C]] <- "blue"
  colVec[scopes$hatUC[, index_C]] <- "red"
  plot(scopes$x, y, col = colVec,
       pch = 18, xlab = xlab, ylab = ylab, main = title)
  if(scopes$SCoPEStype == "extraction"){
    lines(scopes$x, C[,1], lty = 2, col = "blue")
    lines(scopes$x, C[,2], lty = 2, col = "red")
  }else{
    lines(scopes$x, C[,1], lty = 1, col = "orchid3")
  }
  lines(scopes$x, scopes$q * scopes$hatsigma * scopes$tN, lty = 2, col = "orchid3")
  lines(scopes$x, -scopes$q * scopes$hatsigma * scopes$tN, lty = 2, col = "orchid3")
}
#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
sim_SCoPES <- function(Msim, N, alpha, C, q_method, model, I = NULL,
                       inclusion = list(L = "inequal", U = "inequal")){
  # Dimension of the tube defining matrix
  dC    = dim(C)
  if(is.null(dC)){
    dC <- c(length(C), 1)
  }
  # Get standard values for the indices of C^\pm if not specified
  if(is.null(I)){
    if(SCoPEStype %in% c( "extraction", "relevance", "lrelevance",
                          "equivalence", "lequivalence")){
      Iminus <- seq(1, dC[2], 2)
      Iplus  <- seq(2, dC[2], 2)
    }else{
      Iminus <- 1:dC[2]
      Iplus  <- 1:dC[2]
    }
    I = list(minus = Iminus, plus = Iplus)
  }else if(is.list(I)){
    if(length(I) == 2){
      Iminus = I$minus
      Iminus = I$plus
    }else{
      stop("I must be a list with two entries indicating which columns of
           C belong to C^- and which to C^+")
    }
  }else{
    stop("I must be a list with two entries indicating which columns of
           C belong to C^- and which to C^+")
  }

  contain <- 0
  hatq    <- rep(NaN, Msim)
  hatmu1C <- list(minus = matrix(FALSE, length(model$x), Msim),
                  plus  = matrix(FALSE, length(model$x), Msim))
  Lcontain <- Ucontain <- matrix(FALSE, length(model$x), Msim)
  detectL  <- array(FALSE, dim = c(dC[1], length(Iminus), Msim))
  detectU  <- array(FALSE, dim = c(dC[1], length(Iplus), Msim))
  NDetect_BH <- NDetect_hommel <- NDetect_sidak <- NDetect  <-  matrix(NaN, 2, Msim)

  for(m in 1:Msim){
    # Get the sample for the simulation run
    Y = SampleFields::SignalPlusNoise(N,
                                      x = model$x,
                                      mu = model$mu,
                                      noise = model$noise,
                                      sigma = model$sigma)
    Y = Y$values
    # Get the residuals
    R = Y - rowMeans(Y)
    q_method$R = R

    # Estimate the mean and sd
    hatmu = rowMeans(Y)
    if(model$truesigma){
      hatsigma = model$sigma(model$x)
    }else{
      hatsigma = apply(Y, 1, sd)
    }

    if(!is.null(q_method$mu1Cest)){
      if(q_method$mu1Cest == "m0"){
        pvals <- apply(Y, 1, function(v) t.test(x = v,
                                                alternative = "two.sided",
                                                conf.level = 1-alpha)$p.value)
        q_method$m0 = min(2*sum(pvals >= 0.5), length(pvals))
      }
    }

    # Get the SCoPES for the method
    res_m <- SCoPES(alpha = alpha, C = C, x = model$x, hatmu = hatmu,
                    hatsigma = hatsigma, tN = 1 / sqrt(N),
                    method = q_method,
                    inclusion = inclusion, mu = model$mu(model$x))

    # Save the useful variables from the simulation
    hatq[m]    <- res_m$q
    if(!is.null(q_method$mu1Cest)){
      if(q_method$mu1Cest == "m0"){
        hatmu1C$minus[, m] <- rep(T, length(model$x))
        hatmu1C$plus[, m]  <- rep(T, length(model$x))
      }else{
        hatmu1C$minus[, m] <- res_m$hatmu1C$minus
        hatmu1C$plus[, m]  <- res_m$hatmu1C$plus
      }
    }else{
      hatmu1C$minus[, m] <- res_m$hatmu1C$minus
      hatmu1C$plus[, m]  <- res_m$hatmu1C$plus
    }
    detectL[,, m]      <- res_m$hatLC
    detectU[,, m]      <- res_m$hatUC
    Lcontain[, m]      <- res_m$Lcontain_loc
    Ucontain[, m]      <- res_m$Ucontain_loc
    NDetect[1, m]      <- res_m$NtrueDetect
    NDetect[2, m]      <- res_m$NfalseDetect
    contain <- contain + res_m$contain / Msim

    if(SCoPEStype == "classical"){
      I0 = model$mu(model$x) == C
      I1 = model$mu(model$x) != C
      #---------------------------------------------------------------------------
      # Testing alternatives
      pvals <- apply(Y, 1, function(v) t.test(x = v,
                                              alternative = "two.sided",
                                              conf.level = 1-alpha)$p.value)

      # pvals_adjust = p.adjust(pvals, method = "sidak")
      detect_hommel  = MPT_Detect(alpha, pvals, "hommel")
      detect_sidak = FWE_control(alpha, pvals, "sidak")
      detect_BH    = FDR_control(alpha, pvals, "BH")

      NDetect_hommel[1, m]  <- sum(detect_hommel[I1])
      NDetect_hommel[2, m]  <- sum(detect_hommel[I0])
      NDetect_sidak[1, m] <- sum(detect_sidak[I1])
      NDetect_sidak[2, m] <- sum(detect_sidak[I0])
      NDetect_BH[1, m] <- sum(detect_BH[I1])
      NDetect_BH[2, m] <- sum(detect_BH[I0])
      }
  }

    if(SCoPEStype == "classical"){
      return(list(coverage  = contain, Lcoverage = Lcontain,
                  Ucoverage = Ucontain, q = hatq, mu1C = hatmu1C,
                  detectL   = detectL, detectU = detectU,
                  NDetect       = NDetect,
                  NDetect_hommel  = NDetect_hommel,
                  NDetect_sidak = NDetect_sidak,
                  NDetect_BH    = NDetect_BH))
    }else{
      return(list(coverage = contain, Lcoverage = Lcontain,
                  Ucoverage   = Ucontain, q = hatq, mu1C = hatmu1C,
                  detectL     = detectL, detectU = detectU, NDetect = NDetect))
    }
}
