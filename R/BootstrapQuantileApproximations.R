#------------------------------------------------------------------------------#
#                                                                              #
#     Bootstrap functions to obtain the quantiles of the maximum of random     #
#     fields                                                                   #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#      - matrixStats
#
# Contained functions:
#      - ParametricBootstrap (tested)
#      - MultiplierBootstrap (tested)
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline included
#------------------------------------------------------------------------------#

#' Parametric bootstrap estimator for the quantile of the maximum of a random field.
#' The function contains a bootstrap-t version based on Degras (2011)
#' "Simultaneous confidence bands for non-parametric regression with functional data" and
#' a simple non-parametric bootstrap. For large sample sizes the versions do agree, however,
#' for small sample sizes the bootstrap-t has better covering rates, but higher variability
#' in the estimate of the quantile than the simple non-paramtric bootstrap.
#'
#' @param A array of dimension K_1 x ... x K_D x N containing N-realizations of
#'  a random field over a D-dimensional domain.
#' @param alpha numeric the targeted upper quantile of the maxiumum of the
#'   absolute value of the random field. Default is 0.95.
#' @param Mboots numeric amount of bootstrap replicates. Default is 5e3.
#' @param method string specifies the bootstrap version. Options are "t" and
#'  "regular". Default is "t".
#' @return list with elements
#'  \itemize{
#'   \item z Vector containing the realizations of the bootstrap of the approximation of the maximum of the random field
#'   \item q alpha-quantile value of z
#' }
#' @export
ParametricBootstrap <- function( A,
                                 alpha  = 0.05,
                                 Mboots = 5e3,
                                 method = "t" ){

  #---- Check user input and put default values
  # Check A
  if( !is.array( A ) ){
    stop( "'A' must be an array." )
  }

  #---- Check params and put defaults, if missing
  # Check Mboots
    if( is.numeric( Mboots ) ){
      if( Mboots %% 1 != 0 & Mboots <= 0 ){
        stop( "The input 'Mboots' needs to be a positiv natural number." )
      }
    }else{
      stop("The input 'Mboots' needs to be a positiv natural number.")
    }

  # Check alpha
  if( is.numeric( alpha ) ){
    if( alpha <= 0 | alpha >= 1 ){
      stop("The input 'alpha' needs to be a real number between 0 and 1.")
    }
  }else{
    stop("The input 'alpha' needs to be a real number between 0 and 1.")
  }

  # Check method
  if( is.character( method ) ){
    if( !( method %in% c( "t", "regular" ) ) ){
      stop( "Please, specify a valid choice for 'method'. Options are 't' and 'regular'.")
    }
  }else{
    stop( "Please, specify a valid choice for 'method'. Options are 't' and 'regular'.")
  }

  #----- Precompute useful constants
  # dimension of input
  dimA = dim( A )

  # number of samples
  N    = dimA[length(dimA)]
  D    = length(dimA) - 1

  # Get the bootstrap weights. Note that COLUMNS are realizations of rmultinom!
  counter = rmultinom(Mboots, size = N, prob = rep(1 / N, N))

  #----- Switch by dimension for computational speed
  if( D == 1 ){
    # Bootstrapped and original means
    bootMeans <- A %*% counter / N
    meanA     <- rowMeans(A)

    # Estimate the variance from the sample
    if( method == "regular" ){
      data.sigma <- sqrt( matrixStats::rowVars( A ) )
    }else if( method == "t" ){
      bootSecMoments <- A^2 %*% counter / N
      # We put an abs here to make sure that no NaNs are produced due to machine precision error.
      data.sigma <- sqrt( ( N / ( N - 1 ) ) * abs( bootSecMoments - bootMeans^2 ) )
    }

    # Bootstrap realisations of the maximum of the t-field approximation
    distVec <- sqrt(N) * apply(abs((bootMeans - meanA) / data.sigma),
                               2, max )

  }else if( D > 1 ){
    # Bootstrapped and original means
    meanA     <- apply( A, 1:D, mean )
    bootMeans <- array( matrix( A, prod( dimA[ -( D + 1 ) ] ), N ) %*% counter,
                        dim = c( dimA[ -( D + 1 ) ], Mboots ) ) / N

    # Estimate the variance from the sample
      if( method == "regular" ){
        data.sigma <- array( rep( sqrt( apply( A, 1:D, var ) ), Mboots ),
                             dim = c( dimA[ -( D + 1 ) ], Mboots ) )
      }else if( method == "t" ){
        bootSecMoments <- array( matrix( A^2, prod( dimA[ -( D + 1 ) ] ), N ) %*% counter,
                                 dim = c( dimA[ -( D + 1 ) ], Mboots ) ) / N
        # We put an abs here to make sure that no NaNs are produced due to machine precision error.
        data.sigma <- sqrt( ( N / ( N - 1 ) ) * abs( bootSecMoments - bootMeans^2 ) )
      }

    # Compute bootstrap distribution of the absolute value of the maximum
    distVec <- sqrt( N ) * apply(
                abs( ( bootMeans - array( rep( meanA, Mboots ),
                                          dim = c( dimA[ -( D + 1 ) ], Mboots ) ) ) / data.sigma ), D + 1, max)

  }

  # Return list containing the bootstrap values and the estimated quantile
    return( list( z = distVec,
                  q = quantile( distVec, 1 - alpha, type = 8 ) ) )
}

#' Multiplier bootstrap estimator for the quantile of the maximum of a random field. The
#' function contains a multiplier-t version based on Telschow Schwartzmann (2019)
#' "Simultaneous confidence bands for non-parametric regression with functional data"
#' and a regular version based on Cherno. For large sample sizes the versions do agree,
#' however, for small sample sizes the bootstrap-t has better covering rates, but a slightly
#' higher variability in the estimate of the quantile than the simple version based on Cherno.
#'
#' @param R array of dimension K_1 x ... x K_D x N containing N-realizations of
#'  a random field over a D-dimensional domain.
#' @param Q array of dimension K_1 x ... x K_D x N containing N-realizations of
#'  a random field over a D-dimensional domain. Default NULL, i.e. one sample
#'  case.
#' @param alpha numeric the targeted upper quantile of the maxiumum of the
#'   absolute value of the random field. Default is 0.95.
#' @param Mboots numeric amount of bootstrap replicates. Default is 5e3.
#' @param method string specifies the bootstrap version. Options are "t" and
#'  "regular". Default is "t".
#' @param weights string specifying the multipliers. Options are "gaussian",
#'  "rademacher" and "mammen". Default is "rademacher".
#' @return list with elements
#'  \itemize{
#'   \item z Vector containing the realizations of the bootstrap of the approximation of the maximum of the random field
#'   \item q alpha-quantile value of z
#' }
#' @export
MultiplierBootstrap <- function( R,
                                 Q       = NULL,
                                 alpha   = 0.05,
                                 Mboots  = 5e3,
                                 method  = "t",
                                 weights = "rademacher"){

  #---- Check user input and put default values
  # Check R
  if( !is.array( R ) ){
    stop("'R' must be an array.")
  }
  # Check Q
  if( !( is.array( Q ) | is.null( Q ) ) ){
    stop("'Q' must be an array.")
  }

  #---- Check params and put defaults, if missing
  # Check Mboots
  if( is.numeric( Mboots ) ){
    if( Mboots %% 1 != 0 & Mboots <= 0 ){
      stop( "The input 'Mboots' needs to be a positiv natural number." )
    }
  }else{
    stop("The input 'Mboots' needs to be a positiv natural number.")
  }

  # Check alpha
  if( is.numeric( alpha ) ){
    if( alpha <= 0 | alpha >= 1 ){
      stop("The input 'alpha' needs to be a real number between 0 and 1.")
    }
  }else{
    stop("The input 'alpha' needs to be a real number between 0 and 1.")
  }

  # Check method
  if( is.character( method ) ){
    if( !( method %in% c( "t", "regular" ) ) ){
      stop( "Please, specify a valid choice for 'method'. Options are 't' and 'regular'.")
    }
  }else{
    stop( "Please, specify a valid choice for 'method'. Options are 't' and 'regular'.")
  }

  # Check weights
  if( is.character( weights ) ){
    if( !( weights %in% c( "gaussian", "rademacher", "mammen" ) ) ){
      stop( "Please, specify a valid choice for 'weights'. Options are 'gaussian', 'rademacher' and 'mammen'.")
    }
  }else{
    stop( "Please, specify a valid choice for 'weights'. Options are 'gaussian', 'rademacher' and 'mammen'.")
  }

  #----- One sample case
  if( is.null( Q ) ){
      #----- Precompute useful constants
      # dimension of input
      dimR = dim( R )

      # number of samples
      N    = dimR[ length( dimR ) ]
      D    = length( dimR ) - 1

      #----- Simulate multiplier weights
      if( weights == "gaussian" ){
        multiplier <- matrix( rnorm( N * Mboots ), N, Mboots )
      }else if( weights == "rademacher" ){
        multiplier <- matrix( sample( c( -1, 1 ), N * Mboots, replace = T ), N, Mboots )
      }else{
        multiplier <- matrix( sqrt( 5 ) *
                              rbinom( N * Mboots,
                                      1,
                                      ( sqrt( 5 ) - 1 ) / 2 / sqrt( 5 ) ) +
                              ( 1 - sqrt( 5 ) ) / 2,
                              N,
                              Mboots )
      }

      #----- Switch by dimension for computational speed
      if( D == 1 ){
        # Compute bootstrap means
        bootMeans <- R %*% multiplier / N

        # Estimate the variance from the sample
          if( method == "regular" ){
            data.sigma <- sqrt( matrixStats::rowVars( R ) )
          }else if( method == "t" ){
            bootSecMoments <- R^2 %*% multiplier^2 / N
            # We put an abs here to make sure that no NaNs are produced due to machine precision error.
            data.sigma <- sqrt( ( N / ( N - 1 ) ) * abs( bootSecMoments - bootMeans^2 ) )
          }

        # Compute bootstrap distribution of the maximum
        distVec <- sqrt( N ) * apply( abs( bootMeans / data.sigma ), 2, max )

      }else if( D > 1 ){
        # Compute bootstrap means
        bootMeans = array( matrix( R, prod( dimR[ -( D + 1 ) ] ), N ) %*% multiplier,
                                   dim = c( dimR[ -( D + 1 ) ], Mboots ) ) / N

        # Estimate the variance from the sample
        if( method == "regular" ){
          data.sigma <- sqrt( apply( R, 1:D, var  ) )
          distVec    <- apply( array( as.vector( abs( bootMeans ) ) / as.vector( data.sigma ),
                                      dim = c( dimR[ -( D + 1 ) ], Mboots ) ), D + 1, max)  * sqrt( N )
        }else if( method == "t" ){
          bootSecMoments <- array( matrix( R^2, prod( dimR[ -( D + 1 ) ] ), N ) %*% multiplier^2,
                                   dim = c( dimR[ -( D + 1 ) ], Mboots ) ) / N
          # We put an abs here to make sure that no NaNs are produced due to machine precision error.
          data.sigma <- sqrt( ( N / ( N - 1 ) ) * abs( bootSecMoments - bootMeans^2 ) )
          distVec    <- apply( abs( bootMeans ) / data.sigma, D + 1, max ) * sqrt( N )
        }
      }
  }else{
    #----- Precompute useful constants
    # dimension of input
    dimR = dim( R )
    dimQ = dim( Q )
    # number of samples
    N    = dimR[ length( dimR ) ]
    M    = dimQ[ length( dimQ ) ]
    c    = N / M
    D    = length( dimR ) - 1

    #----- Obtain multiplier weights
    if( weights == "gaussian" ){
      multiplierR <- matrix( rnorm( N * Mboots ), N, Mboots )
      multiplierQ <- matrix( rnorm( M * Mboots ), M, Mboots )
    }else{
      multiplierR <- matrix( sample( 1:2, N * Mboots, replace = T ) * 2 - 3, N, Mboots )
      multiplierQ <- matrix( sample( 1:2, M * Mboots, replace = T ) * 2 - 3, M, Mboots )
    }

    #----- Switch by dimension for computational speed
    if( D == 1 ){
      #----- Compute bootstrap means
      bootMeansR <- R %*% multiplierR / N
      bootMeansQ <- Q %*% multiplierQ / M

      #----- Estimate the variance from the samples
      if( method == "regular" ){
        data.sigmaR <- matrixStats::rowVars( R )
        data.sigmaQ <- matrixStats::rowVars( Q )
      }else if( method == "t" ){
        bootSecMoments <- R^2 %*% multiplierR^2 / N
        data.varR      <- ( N / ( N - 1 ) ) * abs( bootSecMoments - bootMeansR^2 )

        bootSecMoments <- Q^2 %*% multiplierQ^2 / M
        data.varQ     <- ( M / ( M - 1 ) ) * abs( bootSecMoments - bootMeansQ^2 )
      }
      data.sigma  <- sqrt( ( 1 + 1 / c ) * data.varR + ( 1 + c ) * data.varQ )

      #----- Compute bootstrap distribution of the maximum
      distVec <- sqrt( N + M ) * apply( abs( ( bootMeansR + bootMeansQ ) / data.sigma ), 2, max )

    }else if( D > 1 ){
      #----- Compute bootstrap means
      bootMeansR = array( matrix( R, prod( dimR[ -( D+1 ) ] ), N ) %*%
                      multiplierR, dim = c( dimR[ -( D+1 ) ], Mboots ) ) / N
      bootMeansQ = array( matrix( Q, prod( dimR[ -( D+1 ) ] ), M ) %*%
                      multiplierQ, dim = c( dimQ[ -( D+1 ) ], Mboots ) ) / M

      #----- Estimate the variance from the sample
      if( method == "regular" ){
        data.varR <- apply( R, 1:D, var )
        data.varQ <- apply( Q, 1:D, var )
        data.sigma  <- sqrt( ( 1 + 1 / c ) * data.varR + ( 1 + c ) * data.varQ )

        distVec    <- apply( array( as.vector( abs( bootMeans ) ) /
                                      as.vector( data.sigma ),
                                    dim = c( dimR[ -( D+1 ) ], Mboots ) ), D + 1, max ) * sqrt( N + M - 2 )
      }else if( method == "t" ){
        #----- bootstrapped variance for R
        bootSecMoments <- array( matrix( R^2, prod( dimR[ -( D+1 ) ] ), N ) %*%
                                   multiplierR^2,
                                 dim = c( dimR[ -( D+1 ) ], Mboots ) ) / N
        data.varR <- ( N / ( N - 1 ) ) * abs( bootSecMoments - bootMeansR^2 )

        #----- bootstrapped variance for Q
        bootSecMoments <- array( matrix( Q^2, prod( dimQ[ -( D+1 ) ] ), N ) %*%
                                   multiplierQ^2,
                                 dim = c( dimQ[ -( D+1 ) ], Mboots ) ) / M
        data.varQ <- ( M / ( M - 1 ) ) * abs( bootSecMoments - bootMeansQ^2 )

        #----- bootstrapped variance of the asymptotic process
        data.sigma  <- sqrt( ( 1 + 1 / c ) * data.varR + ( 1 + c ) * data.varQ )

        #----- max statistic
        distVec <- apply( abs( bootMeansR + bootMeansQ )  / data.sigma,
                          D+1,
                          max )  * sqrt(N+M)
      }
    }
  }

  q = quantile( distVec, 1 - alpha, type = 8 )

  #----- Return quantile and bootstrap distribution
  return( list( z       = distVec,
                q       = q,
                samples = bootMeans / data.sigma * sqrt( N ) ) )
}


#' Multiplier bootstrap estimator for the quantile of the maximum of a random field. The
#' function contains a multiplier-t version based on Telschow Schwartzmann (2019)
#' "Simultaneous confidence bands for non-parametric regression with functional data"
#' and a regular version based on Cherno. For large sample sizes the versions do agree,
#' however, for small sample sizes the bootstrap-t has better covering rates, but a slightly
#' higher variability in the estimate of the quantile than the simple version based on Cherno.
#'
#' @param R array of dimension K_1 x ... x K_D x N containing N-realizations of
#'  a random field over a D-dimensional domain.
#' @param Q array of dimension K_1 x ... x K_D x N containing N-realizations of
#'  a random field over a D-dimensional domain. Default NULL, i.e. one sample
#'  case.
#' @param alpha numeric the targeted upper quantile of the maxiumum of the
#'   absolute value of the random field. Default is 0.95.
#' @param Mboots numeric amount of bootstrap replicates. Default is 5e3.
#' @param method string specifies the bootstrap version. Options are "t" and
#'  "regular". Default is "t".
#' @param weights string specifying the multipliers. Options are "gaussian",
#'  "rademacher" and "mammen". Default is "rademacher".
#' @return list with elements
#'  \itemize{
#'   \item z Vector containing the realizations of the bootstrap of the approximation of the maximum of the random field
#'   \item q alpha-quantile value of z
#' }
#' @export
MultiplierBootstrapSplit <- function( R,
                                 Splus,
                                 Sminus,
                                 Q = NULL,
                                 alpha   = 0.05,
                                 Mboots  = 5e3,
                                 method  = "t",
                                 weights = "rademacher" ){

  #---- Check user input and put default values
  # Check R
  if( !is.array( R ) ){
    stop("'R' must be an array.")
  }
  # Check Q
  if( !( is.array( Q ) | is.null( Q ) ) ){
    stop("'Q' must be an array.")
  }

  #---- Check params and put defaults, if missing
  # Check Mboots
  if( is.numeric( Mboots ) ){
    if( Mboots %% 1 != 0 & Mboots <= 0 ){
      stop( "The input 'Mboots' needs to be a positiv natural number." )
    }
  }else{
    stop("The input 'Mboots' needs to be a positiv natural number.")
  }

  # Check alpha
  if( is.numeric( alpha ) ){
    if( alpha <= 0 | alpha >= 1 ){
      stop("The input 'alpha' needs to be a real number between 0 and 1.")
    }
  }else{
    stop("The input 'alpha' needs to be a real number between 0 and 1.")
  }

  # Check method
  if( is.character( method ) ){
    if( !( method %in% c( "t", "regular" ) ) ){
      stop( "Please, specify a valid choice for 'method'. Options are 't' and 'regular'.")
    }
  }else{
    stop( "Please, specify a valid choice for 'method'. Options are 't' and 'regular'.")
  }

  # Check weights
  if( is.character( weights ) ){
    if( !( weights %in% c( "gaussian", "rademacher", "mammen" ) ) ){
      stop( "Please, specify a valid choice for 'weights'. Options are 'gaussian', 'rademacher' and 'mammen'.")
    }
  }else{
    stop( "Please, specify a valid choice for 'weights'. Options are 'gaussian', 'rademacher' and 'mammen'.")
  }

  #----- One sample case
    #----- Precompute useful constants
    # dimension of input
    dimR = dim( R )

    # number of samples
    N    = dimR[ length( dimR ) ]
    D    = length( dimR ) - 1

    #----- Simulate multiplier weights
    if( weights == "gaussian" ){
      multiplier <- matrix( rnorm( N * Mboots ), N, Mboots )
    }else if( weights == "rademacher" ){
      multiplier <- matrix( sample( c( -1, 1 ), N * Mboots, replace = T ), N, Mboots )
    }else{
      multiplier <- matrix( sqrt( 5 ) *
                              rbinom( N * Mboots,
                                      1,
                                      ( sqrt( 5 ) - 1 ) / 2 / sqrt( 5 ) ) +
                              ( 1 - sqrt( 5 ) ) / 2,
                            N,
                            Mboots )
    }

      # Compute bootstrap means
      bootMeans <- R %*% multiplier / N

      # Estimate the variance from the sample
      if( method == "regular" ){
        data.sigma <- t(t(sqrt( matrixStats::rowVars( R ) )))
      }else if( method == "t" ){
        bootSecMoments <- R^2 %*% multiplier^2 / N
        # We put an abs here to make sure that no NaNs are produced due to machine precision error.
        data.sigma <- sqrt( ( N / ( N - 1 ) ) * abs( bootSecMoments - bootMeans^2 ) )
      }

      # Compute bootstrap distribution of the maximum
      if(sum(Splus) == 1){
        distVec1 <- sqrt( N ) * bootMeans[Splus,] / data.sigma[Splus,]
      }else if(sum(Splus) > 1){
        distVec1 <- sqrt( N ) * apply( t(t(bootMeans[Splus,])) / data.sigma[Splus,], 2, max )
      }else{
        distVec1 <- rep(-Inf, Mboots)
      }
      if(sum(Sminus) == 1){
        distVec2 <- -sqrt( N ) * bootMeans[Sminus,] / data.sigma[Sminus,]
      }else if(sum(Sminus) > 1){
        distVec2 <- -sqrt( N ) * apply( t(t(bootMeans[Sminus,])) / data.sigma[Sminus,], 2, max )
      }else{
        distVec2 <- rep(-Inf, Mboots)
      }

      if( method == "regular" ){
        data.sigma = as.vector(data.sigma)
      }

      distVec  <- apply( rbind(distVec1, distVec2), 2, max )

  #----- Return quantile and bootstrap distribution
  return( list( z       = distVec,
                q       = quantile( distVec, 1 - alpha, type = 8 ),
                samples = bootMeans / data.sigma * sqrt( N ) ) )
}
