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
#' @param coords
#' @param mask
#' @return list with elements
#'  \itemize{
#'   \item hatmean pointwise sample mean
#'   \item scb list containing the upper and lower bounds of the simultaneous confidence band
#'   \item level targeted covering probability
#'   \item q quantile of the maximum of the residual field
#' }
#' @export
scb_moments <- function(  Y,
                          level          = .95,
                          transformation = "linear",
                          moments        = NULL,
                          method  = list(name = "tGKF", field = "t"),
                          bias    = FALSE,
                          se.est  = TRUE,
                          coords  = NULL,
                          mask    = NULL ){
  #----- Check user input
  # Check input Y
  if( !is.array( Y ) ){
    stop( "Y must be an array containing the realisations of
           the functional data, i.e. an array with last
           dimension enumerating the realisations." )
  }

  # Check input level
  if( is.numeric( level ) ){
    if( level >= 1 | level <= 0 ){
      stop( "The input 'level' needs to be strictly between 0 and 1." )
    }
  }else{
    stop( "The input 'level' needs to be strictly between 0 and 1." )
  }

  # Check input method
  if( is.character( method$name ) ){
    if( !( method$name %in% c( "tGKF", "GKF", "ffscb",
                               "ParamBoot", "MultBoot") ) ){
      stop( "Choose a valid option from the available quantile approximations.
            Please, check the help page." )
    }
  }else{
      stop( "Choose a valid option from the available quantile approximations.
            Please, check the help page." )
  }

  #----- Compute constants
  dimY = dim( Y )
  # dimension of domain
  D = length( dimY ) - 1
  # get number of sample curves
  N = dimY[ length( dimY ) ]

  #----- Get the correct residuals
  residuals = DeltaMomentResiduals( Y = Y,
                                    transformation = transformation,
                                    moments = moments )

  # normalize the residuals, if neccessary
  if( method$name == "ParamBoot" ){
    normResiduals = Y
  }else{
    # Get the normalized residuals
    normResiduals = residuals$delta.res /
                                  array( residuals$delta.sd, dim = c( dimY ) )
  }

  #----- Compute the maximum quantile
  q = maxQuantile( R      = normResiduals,
                   alpha  = 1 - level,
                   method = method,
                   coords = coords,
                   mask   = mask )

  #----- Compute the SCB
  # Get the se estimate or use the true for Gaussian data
  if( se.est ){
    se  <- sqrt( N / ( N - 1 ) ) * residuals$delta.sd / sqrt( N )
  }else{
    if( transformation == "cohensd" ){
      se <- sqrt( 1 + residuals$statistic^2 / 2 ) / sqrt(N)
    }else if( transformation == "skewness" ){
      se <- sqrt( 6 * ( N-2 ) / ( N+1 ) / ( N+3 ) )
    }else if( transformation %in% c( "kurtosis", "kurtosis (unbiased)" ) ){
      se <- sqrt( 24 * N * (N-2) * (N-3) / (N+1)^2 / (N+3) / (N+5) )
    }
  }

  if( bias ){
    scb <- residuals$statistic - residuals$hatbias - q *  se
    scb <- cbind( scb, residuals$statistic - residuals$hatbias )
    scb <- cbind( scb, residuals$statistic - residuals$hatbias + q * se )
  }else{
    scb <- residuals$statistic - q *  se
    scb <- cbind( scb, residuals$statistic )
    scb <- cbind( scb, residuals$statistic + q * se )
  }

  colnames( scb ) <- c( "SCB.low", "est", "SCB.up" )

  #----- Return a list containing estimate of the mean, SCBs etc.
  if( bias ){
    list( scb   = scb,
          level = level,
          q     = q,
          residuals = residuals,
          bias  = residuals$hatbias  )
  }else{
    list( scb   = scb,
          level = level,
          q     = q,
          residuals = residuals  )
  }
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
maxQuantile <- function( R,
                         alpha,
                         method,
                         coords,
                         mask ){

  #----- different approximation of the max quantile
  if( method$name %in% c( "tGKF", "GKF" ) ){
    #----- Check whether the correct fields are in the list
    # euler characteristic of domain
    if( is.null( method$L0 ) ){
      method$L0 = 1
    }

    # estimator of LKCs
    if( is.null( method$LKC_estim ) ){
      method$LKC_estim <- function( R ){
                                    RFT::LKC_integral( R      = R,
                                                       coords = coords,
                                                       mask   = mask )
                          }
    }

    # marginal distribution of field
    if( is.null( method$field ) ){
      method$field <- "t"
    }

    # degrees of freedom
    if( is.null( method$df ) ){
      method$df <- dim(R)[ length( dim( R ) ) ] - 1
    }

    #----- Get the quantiles of the maximum
    # Estimate the LKCs from the normed residuals
    LKC = c( method$L0,
             method$LKC_estim( R ) )

    # Estimate the quantile using the GKF approximation
    q   = RFT::GKFthreshold( alpha = alpha / 2,
                             LKC   = LKC,
                             type  = method$field,
                             df    = method$df,
                             interval = c( 0, 200 ) )$threshold
  } else if( method$name == "ParamBoot" ){
    # add default value Mboots
    if( is.null( method$Mboots ) ){
      method$Mboots <- 5e4
    }

    # add default value method
    if( is.null( method$method ) ){
      method$method <- "t"
    }

    # add the level to the parameters of the bootstrap
    method$alpha = alpha

    # Estimate the quantile
    q <- ParametricBootstrap( R,
                              alpha  = alpha,
                              Mboots = method$Mboots,
                              method = "t" )$q

  }else if( method$name == "MultBoot" ){
    # add default value Mboots
    if( is.null( method$Mboots ) ){
      method$Mboots <- 5e4
    }

    # add default value method
    if( is.null( method$method ) ){
      method$method <- "t"
    }

    # add default value multiplier
    if( is.null( method$weights ) ){
      method$weights <- "rademacher"
    }

    # estimate the quantile
    q <- MultiplierBootstrap( R,
                              Q       = NULL,
                              alpha   = alpha,
                              Mboots  = method$Mboots,
                              method  = method$method,
                              weights = method$weights )$q
  }

  return( q )
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
                          bias,
                          se.est,
                          trueValue,
                          ... ){
  # Generate output variable
  covRate <- matrix( 0, length( Nvec ), length( method ) )
  rownames( covRate ) <- Nvec
  colnames( covRate ) <- names(method)

  # Monte Carlo simulation loop
  for( nn in 1:length(Nvec) ){
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
        scb = scb_moments( Y              = Y1$values,
                           level          = level,
                           transformation = transformation,
                           moments        = moments,
                           method         = method[[k]],
                           bias           = bias,
                           se.est         = se.est,
                           coords         = list( x = Y1$locations ),
                           mask           = NULL )
        if( all( scb$scb[,1] < trueValue ) & all( trueValue <= scb$scb[,3] ) ){
          covRate[nn,k] = covRate[nn,k] + 1 / Msim
        }
      }
    }
    #print(covRate)
  }
  return( covRate )
}
