#------------------------------------------------------------------------------#
#                                                                              #
#     Functions to compute simultaneous confidence bands for functional data   #
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

#' Computes simultaneous confidence bands for the mean of a sample from a one
#' dimensional functional signal plus noise model. It is possible to choose
#' between different estimators for the quantile.
#'
#' @param Y array of dimension K_1 x ... x K_D x N containing N-realizations of
#' a random field over a D-dimensional domain.
#' @param level numeric the targeted covering probability. Default 0.95.
#' @param transformation expression of the form f ~ f(mu1, ..., muk ) giving
#' the transformation of the moments. Can be a "linear", "cohensd", "skewness" or
#' "kurtosis". Then standard estimators are used.
#' Current options are "skewness" and "kurtosis".
#' @param moments string of the form c("mu1", ..., "muk") containing the names of
#' the moments contained in the order of transformation input.
#' @param method string specifying the method to construct the scb, i.e.
#' estimatate the quantile.
#' Current options are "tGKF", "GKF", "NonparametricBootstrap", "MultiplierBootstrap".
#' Default value is "tGKF".
#' @param bias.est a string indicating the bias estimation method.
#'  \itemize{
#'   \item "asymptotic gaussian" means bias = 0
#'   \item "exact gaussian" uses the true finite sample variance bias assuming
#' the transdormation is cohensd, skewness or kurtosis
#'   \item "estimate" means that the bias is estimated from the delta residuals
#' }
#' Default is bias = 0, i.e., "asymptotic gaussian".
#' ,
#' @param se.est
#'  \itemize{
#'   \item "asymptotic gaussian" means true Gaussian asymptotic variance is used in the
#'   standard error, i.e. se("skewness") = 6 / sqrt(sample size)
#'   \item "exact gaussian" uses the true finite sample variance bias assuming
#' the transdormation is cohensd, skewness or kurtosis
#'   \item "estimate" means that the bias is estimated from the delta residuals
#' }
#' @param smoothing either NULL indicating "no smoothing applied" or a list.
#' In the latter case locpoly from KernSmooth is used
#'  \itemize{
#'   \item If a field "bandwidth" is contained the value is used for the
#'   smoothing. Othererwise crossvalidation estimates the optimal bandwidth
#'   \item  If a field "degree" is contained the value is used for the
#'   smoothing. Othererwise degree = 1 is used.
#'    }
#' @param coords test
#' @param mask test
#' @return list with elements
#'  \itemize{
#'   \item hatmean pointwise sample mean
#'   \item scb list containing the upper and lower bounds of the simultaneous confidence band
#'   \item level targeted covering probability
#'   \item q quantile of the maximum of the residual field
#' }
#' @export
scb_moments <- function(Y,
                        level          = .95,
                        transformation = "linear",
                        moments        = NULL,
                        method    = list(name = "GKF", field = "t"),
                        bias.est  = "asymptotic gaussian",
                        se.est    = "estimate",
                        smoothing = NULL,
                        coords    = NULL,
                        mask      = NULL){
  #----- Check user input
  # # Check input Y
  # if( !is.array( Y ) ){
  #   stop( "Y must be an array containing the realisations of
  #          the functional data, i.e. an array with last
  #          dimension enumerating the realisations." )
  # }

  # Check input level
  if(is.numeric(level)){
    if(level >= 1 | level <= 0){
      stop( "The input 'level' needs to be strictly between 0 and 1." )
    }
  }else{
    stop( "The input 'level' needs to be strictly between 0 and 1." )
  }

  # Check input method
  if(is.character(method$name)){
    if(!(method$name %in% c("tGKF", "GKF", "ffscb", "ParamBoot", "MultBoot"))){
      stop( "Choose a valid option from the available quantile approximations.
            Please, check the help page." )
    }
  }else{
      stop( "Choose a valid option from the available quantile approximations.
            Please, check the help page." )
  }
  x = Y$locations
  Y = Y$values
  #----- Compute constants
  dimY = dim(Y)
  # dimension of domain
  D = length(dimY) - 1
  # get number of sample curves
  N = dimY[length(dimY)]

  # define smooth weights
  SmoothWeights <- NULL

  #----- Smooth the data
  if(!is.null(smoothing)){
    if(is.null(smoothing$degree)){
      smoothing$degree = 1
    }
    if(is.null(smoothing$bandwidth)){
      smoothing$bandwidth = SCBmeanfd::cv.select(x = x, y = t(Y), degree = smoothing$degree)
      #smoothing$bandwidth = SCBmeanfd::plugin.select(x = x, y = t(Y), degree = smoothing$degree)
    }
    if(is.null(smoothing$xeval)){
      dx = range(x)
      smoothing$xeval = seq(dx[1], dx[2], length.out = diff(dx) * 175 )
    }
    if(is.null(smoothing$kernel)){
      smoothing$kernel = locpol::gaussK
    }
    # Get the local poly smoothing weights
    SmoothWeights <- as.matrix(
                      locpol::locLinWeightsC(x      = x,
                                             xeval  = smoothing$xeval,
                                             bw     = smoothing$bandwidth,
                                             kernel = smoothing$kernel )$locWeig )

    #### Compute smoothed data
    Y <- SmoothWeights %*% Y
    # put the xeval the new x
    x <- smoothing$xeval
    # change the dimY
    dimY = dim(Y)
  }

  #----- Get the correct residuals
  residuals = DeltaMomentResiduals(Y = Y,
                                   transformation = transformation,
                                   moments = moments)

  # normalize the residuals, if neccessary
  if( method$name == "ParamBoot" ){
    normResiduals = Y
  }else{
    # Get the normalized residuals
    normResiduals = residuals$delta.res /
                                  array(residuals$delta.sd, dim = c(dimY))
  }

  #----- Compute the maximum quantile
  q = maxQuantile( R      = normResiduals,
                   alpha  = 1 - level,
                   method = method,
                   coords = coords,
                   mask   = mask )

  #----- Compute the SCB
  # Get the se estimate or use the true for Gaussian data
  if(se.est == "estimate"){
    se  <- sqrt(N / (N - 1)) * residuals$delta.sd / sqrt(N)
  }else if(se.est == "asymptotic gaussian"){
    if(transformation == "cohensd"){
      se <- sqrt(1 + residuals$statistic^2 / 2) / sqrt(N)
    }else if(transformation == "skewness"){
      se <- 6 / sqrt(N)
    }else if(transformation %in% c("kurtosis", "kurtosis (unbiased)")){
      se <- 24 / sqrt(N)
    }
  }else if(se.est == "exact gaussian"){
      se <- se_Gaussian(transformation, N, hatd = residuals$statistic)
  }else{
    se <- se.est / sqrt(N)
  }

  # Get the se estimate or use the true for Gaussian data
  if(bias.est == "estimate"){
    bias  <- residuals$hatbias
  }else if(bias.est == "asymptotic gaussian"){
    bias <- 0
  }else if(bias.est == "exact gaussian"){
    bias <- bias_Gaussian(transformation, N, hatd = residuals$statistic)
  }

  # Get the confidence bands
  scb <- residuals$statistic - bias - q *  se
  scb <- cbind(scb, residuals$statistic - bias)
  scb <- cbind(scb, residuals$statistic - bias + q * se)

  colnames( scb ) <- c("SCB.low", "est", "SCB.up")

  #----- Return a list containing estimate of the mean, SCBs etc.
  list(scb   = scb,
       level = level,
       q     = q,
       residuals = residuals,
       bias  = bias,
       se    = se,
       Y     = Y,
       locations     = x,
       SmoothWeights = SmoothWeights)
}

#' Computes the quantile of the maximum statistic
#'
#' @param R array of dimension K_1 x ... x K_D x N containing N-realizations of
#' a random field over a D-dimensional domain.
#' @param alpha numeric the upper tail probability
#' @param method
#'
#' @return quantile of the maximum process
#' @export
maxQuantile <- function(R,
                        alpha,
                        method,
                        coords,
                        mask){

  #----- different approximation of the max quantile
  if(method$name == "GKF"){
    #----- Check whether the correct fields are in the list
    # euler characteristic of domain
    if(is.null(method$L0)){
      method$L0 = 1
    }

    # estimator of LKCs
    if(is.null(method$LKC_estim)){
      method$LKC_estim <- function(R){
                                    RFT::LKC_integral(R      = R,
                                                      coords = coords,
                                                      mask   = mask )
                          }
    }

    # marginal distribution of field
    if(is.null(method$field)){
      method$field <- "t"
    }

    # degrees of freedom
    if(is.null(method$df)){
      method$df <- dim(R)[length(dim(R))] - 1
    }

    #----- Get the quantiles of the maximum
    # Estimate the LKCs from the normed residuals
    LKC = c(method$L0, method$LKC_estim(R))

    # Estimate the quantile using the GKF approximation
    q   = RFT::GKFthreshold( alpha = alpha / 2,
                             LKC   = LKC,
                             type  = method$field,
                             df    = method$df,
                             interval = c(0, 200))$threshold
  }else if(method$name == "ParamBoot"){
    # add default value Mboots
    if( is.null( method$Mboots ) ){
      method$Mboots <- 5e4
    }

    # add default value method
    if(is.null(method$method)){
      method$method <- "t"
    }

    # add the level to the parameters of the bootstrap
    method$alpha = alpha

    # Estimate the quantile
    q <- ParametricBootstrap(R,
                             alpha  = alpha,
                             Mboots = method$Mboots,
                             method = "t" )$q

  }else if(method$name == "MultBoot"){
    # add default value Mboots
    if(is.null(method$Mboots)){
      method$Mboots <- 5e4
    }

    # add default value method
    if(is.null(method$method)){
      method$method <- "t"
    }

    # add default value multiplier
    if(is.null(method$weights)){
      method$weights <- "rademacher"
    }

    # estimate the quantile
    q <- MultiplierBootstrap(R,
                             Q       = NULL,
                             alpha   = alpha,
                             Mboots  = method$Mboots,
                             method  = method$method,
                             weights = method$weights )$q
  }

  return(q)
}

#' Computes the quantile of the maximum statistic
#'
#' @param Msim number of simulations.
#' @param Nvec vector containing the sample sizes
#' @param x
#' @param mu
#' @param sigma
#' @param noise
#' @param trueValue
#' @inheritParams scb_moments
#'
#' @return quantile of the maximum process
#' @export
#'
covering_scb <- function( Msim  = 5e4,
                          Nvec,
                          level,
                          method = list( name = "GKF" ),
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
                          ... ){
  # Generate output variable
  covRate <- matrix(0, length(Nvec), length(method))
  rownames(covRate) <- Nvec
  colnames(covRate) <- names(method)

  if(is.function(trueValue)){
    trueValuef <- trueValue
  }else{
    trueValuef <- NULL
  }

  if(suppressWarnings(!is.na(as.numeric(se.est)))){
    se.est = as.numeric(se.est)
  }

  # Monte Carlo simulation loop
  for( nn in 1:length(Nvec) ){
    if(is.function(trueValuef)){
      trueValue <- trueValuef(x, N)
    }
    for( m in 1:Msim ){
      # Generate the data
      Y1 <- SignalPlusNoise( N  = Nvec[nn],
                             x  = x,
                             mu = mu,
                             noise = noise,
                             sigma = sigma, ... )

      for( k in 1:length(method) ){
        # construct the
        N  = Nvec[nn]
        scb = scb_moments( Y              = Y1,
                           level          = level,
                           transformation = transformation,
                           moments        = moments,
                           method         = method[[k]],
                           bias.est       = bias.est,
                           se.est         = se.est,
                           smoothing      = smoothing,
                           coords         = list(x = Y1$locations),
                           mask           = NULL )

        if( all(scb$scb[,1] < trueValue) & all(trueValue <= scb$scb[,3]) ){
          covRate[nn,k] = covRate[nn,k] + 1 / Msim
        }
      }
    }
   # print(covRate)
  }
  return(list(rates = covRate, level = level))
}
