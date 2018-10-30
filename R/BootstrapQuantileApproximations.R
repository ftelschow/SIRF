################################################################################################
##                                                                                            ##
##  This file contains estimators based on the bootstrap to estimate the quantile of the      ##
##  maximum of random fields                                                                  ##
##                                                                                            ##
################################################################################################
## required packages:
##
##
##  included functions:
##     - NonParametricBootstrap (tested)
##     - MultiplierBootstrap (tested)
##
################################################################################################
#' Non-parametric bootstrap estimator for the quantile of the maximum of a random field. The function contains a bootstrap-t version based on Degras (2011) "Simultaneous confidence bands for non-parametric regression with functional data" and a simple non-parametric bootstrap. For large sample sizes the versions do agree, however, for small sample sizes the bootstrap-t has better covering rates, but higher variability in the estimate of the quantile than the simple non-paramtric bootstrap.
#'
#' @param A Array the last dimension enumerates the realizations of the random field.
#' @param params list containing elements
#' \itemize{
#'   \item Mboots Numeric the amount of bootstrap replicates. Default is 5e3.
#'   \item alpha Numeric the targeted upper quantile of the maxiumum of the absolute value of the random field. Default is 0.95.
#'   \item method String the bootstrap version. Currently the options are "t" for the bootstrap-t and "regular". Default is "t".
#'   \item stat String specifying the statistic. Currently the options "mean" (default) and "SNR" are supported.
#' }
#' @return list with elements
#'  \itemize{
#'   \item z Vector containing the realizations of the bootstrap of the approximation of the maximum of the random field
#'   \item q alpha-quantile value of z
#' }
#' @export
NonParametricBootstrap <- function(A, params=NULL){
  ###### Check user input and put default values
  #### Check A
  if( !is.array(A) ){
    stop("'A' must be an array.")
  }
  #### Check params and put defaults, if missing
  ## Check Mboots
  if( is.null(params$Mboots) ){
      Mboots = 5e3
  }else{
    if( is.numeric(params$Mboots) ){
      if( params$Mboots%%1 !=0 & params$Mboots<=0 ){
        stop("The input 'params$Mboots' needs to be a positiv natural number.")
      }else{
        Mboots = params$Mboots
      }
    }else{
      stop("The input 'params$Mboots' needs to be a positiv natural number.")
    }
  }

  ## Check params$alpha
  if( is.null(params$alpha) ){
    alpha = 0.95
  }else{
    if( is.numeric(params$alpha) ){
      if( params$alpha <= 0 | params$alpha >= 1 ){
        stop("The input 'params$alpha' needs to be a real number between 0 and 1.")
      }else{
        alpha = params$alpha
      }
    }else{
      stop("The input 'params$alpha' needs to be a real number between 0 and 1.")
    }
  }
  ## Check params$method
  if( is.null(params$method) ){
    method = "t"
  }else{
    if( is.character(params$method) ){
      if( !(params$method %in% c("t", "regular")) ){
        stop("Please, specify a valid choice for 'params$method'. Options are 't' and 'regular'.")
      }else{
        method = params$method
      }
    }else{
      stop("Please, specify a valid choice for 'params$method'. Options are 't' and 'regular'.")
    }
  }
  if( is.null(params$stat) ){
    stat = "mean"
  }
  ###### Precompute useful constants
  ## dimension of input
  dimA = dim(A)
  ## number of samples
  N    = dimA[length(dimA)]
  D    = length(dimA)-1
  ## Simulate bootstrap weights. Note that COLUMNS are realizations of rmultinom!
  counter = rmultinom( Mboots, size=N, prob=rep(1/N, N) )

  if(stat=="mean"){
      ##### Split the cases into data with 1D domain and higher dimensional domain for computation speed
      if( D==1 ){
        ### Compute means of bootstrap realizations and of the original data.
        bootMeans <- A %*% counter / N
        meanA     <- rowMeans(A)

        ### Estimate the variance from the sample if it is unknown
        if( method == "regular" ){
          data.sigma <- sqrt(matrixStats::rowVars(A))
        }else if( method == "t" ){
          bootSecMoments <- A^2 %*% counter / N
          # We put an abs here to make sure that no NaNs are produced due to machine precision error.
          data.sigma <- sqrt((N / (N-1)) * abs(bootSecMoments - bootMeans^2))
        }

        ### Compute the bootstrap realisation of the maximum of the T-field approximation
        distVec <- sqrt(N) * apply( abs((bootMeans - meanA) / data.sigma), 2, max)

      }else if( D>1 ){
        ### Compute original sample mean and bootstrap means
        meanA <- apply( A, 1:D, mean )
        bootMeans = array( matrix(A, prod(dimA[-(D+1)]), N)%*%counter, dim=c(dimA[-(D+1)], Mboots) ) / N

        ### Estimate the variance from the sample if it is unknown
          if( method == "regular" ){
            data.sigma <- array(rep( sqrt( apply( A, 1:D, var ) ), Mboots ), dim=c(dimA[-(D+1)], Mboots) )
          }else if( method == "t" ){
            bootSecMoments <- array( matrix(A^2, prod(dimA[-(D+1)]), N)%*%counter, dim=c(dimA[-(D+1)], Mboots) ) / N
            # We put an abs here to make sure that no NaNs are produced due to machine precision error.
            data.sigma <- sqrt((N / (N-1)) * abs(bootSecMoments - bootMeans^2))
          }

        ### Compute bootstrap distribution of the absolute value of the maximum
        distVec <- sqrt(N) * apply( abs((bootMeans - array(rep(meanA, Mboots), dim=c(dimA[-(D+1)], Mboots) )) / data.sigma), D+1, max)

      }
  }else if(stat=="SNR"){
    ##### Split the cases into data with 1D domain and higher dimensional domain for computation speed
    if( D==1 ){
      ### Compute means of bootstrap realizations and of the original data.
      bootMeans <- A %*% counter / N
      # We put an abs here to make sure that no NaNs are produced due to machine precision error.
      data.sigma <- sqrt(matrixStats::rowVars(A))

      data.sigma <- sqrt((N / (N-1)) * abs(bootSecMoments - bootMeans^2))


      meanA     <- rowMeans(A)

      ### Estimate the variance from the sample if it is unknown
      if( method == "regular" ){

      }else if( method == "t" ){
        bootSecMoments <- A^2 %*% counter / N

      }

      ### Compute the bootstrap realisation of the maximum of the T-field approximation
      distVec <- sqrt(N) * apply( abs((bootMeans - meanA) / data.sigma), 2, max)

    }else if( D>1 ){
      ### Compute original sample mean and bootstrap means
      meanA <- apply( A, 1:D, mean )
      bootMeans = array( matrix(A, prod(dimA[-(D+1)]), N)%*%counter, dim=c(dimA[-(D+1)], Mboots) ) / N

      ### Estimate the variance from the sample if it is unknown
      if( method == "regular" ){
        data.sigma <- array(rep( sqrt( apply( A, 1:D, var ) ), Mboots ), dim=c(dimA[-(D+1)], Mboots) )
      }else if( method == "t" ){
        bootSecMoments <- array( matrix(A^2, prod(dimA[-(D+1)]), N)%*%counter, dim=c(dimA[-(D+1)], Mboots) ) / N
        # We put an abs here to make sure that no NaNs are produced due to machine precision error.
        data.sigma <- sqrt((N / (N-1)) * abs(bootSecMoments - bootMeans^2))
      }

      ### Compute bootstrap distribution of the absolute value of the maximum
      distVec <- sqrt(N) * apply( abs((bootMeans - array(rep(meanA, Mboots), dim=c(dimA[-(D+1)], Mboots) )) / data.sigma), D+1, max)
    }
  }

  ### Return list containing the bootstrap values and the estimated quantile
  return( list( z = distVec, q = quantile(distVec, 1-alpha, type=8) ) )
}


#' Multiplier bootstrap estimator for the quantile of the maximum of a random field. The function contains a multiplier-t version based on Telschow Schwartzmann (2018) "Simultaneous confidence bands for non-parametric regression with functional data" and a regular version based on Cherno. For large sample sizes the versions do agree, however, for small sample sizes the bootstrap-t has better covering rates, but a slightly higher variability in the estimate of the quantile than the simple version based on Cherno.
#'
#' @param R Array the last dimension enumerates the realizations of an residual field.
#' @param params list containing elements
#' \itemize{
#'   \item Mboots Numeric the amount of bootstrap replicates. Default is 5e3.
#'   \item alpha Numeric the targeted upper quantile of the maxiumum of the absolute value of the random field. Default is 0.95.
#'   \item method String the bootstrap version. Currently the options are "t" for the bootstrap-t and "regular". Default is "t".
#' }
#' @return list with elements
#'  \itemize{
#'   \item z Vector containing the realizations of the bootstrap of the approximation of the maximum of the random field
#'   \item q alpha-quantile value of z
#' }
#' @export
MultiplierBootstrap <- function( R, params=NULL ){
  ###### Check user input and put default values
  #### Check A
  if( !is.array(R) ){
    stop("'A' must be an array.")
  }
  #### Check params and put defaults, if missing
  ## Check Mboots
  if( is.null(params$Mboots) ){
    Mboots = 5e3
  }else{
    if( is.numeric(params$Mboots) ){
      if( params$Mboots%%1 !=0 & params$Mboots<=0 ){
        stop("The input 'params$Mboots' needs to be a positiv natural number.")
      }else{
        Mboots = params$Mboots
      }
    }else{
      stop("The input 'params$Mboots' needs to be a positiv natural number.")
    }
  }

  ## Check params$alpha
  if( is.null(params$alpha) ){
    alpha = 0.95
  }else{
    if( is.numeric(params$alpha) ){
      if( params$alpha <= 0 | params$alpha >= 1 ){
        stop("The input 'params$alpha' needs to be a real number between 0 and 1.")
      }else{
        alpha = params$alpha
      }
    }else{
      stop("The input 'params$alpha' needs to be a real number between 0 and 1.")
    }
  }
  ## Check params$method
  if( is.null(params$method) ){
    method = "t"
  }else{
    if( is.character(params$method) ){
      if( !(params$method %in% c("t", "regular")) ){
        stop("Please, specify a valid choice for 'params$method'. Options are 't' and 'regular'.")
      }else{
        method = params$method
      }
    }else{
      stop("Please, specify a valid choice for 'params$method'. Options are 't' and 'regular'.")
    }
  }

  ###### Precompute useful constants
  # dimension of input
  dimR = dim(R)
  # number of samples
  N    = dimR[length(dimR)]
  D    = length(dimR)-1

  ##### Simulate Gaussian multiplier weights
  multiplier <- matrix(rnorm(N * Mboots), N, Mboots)

  ##### Split the cases into data with 1D domain and higher dimensional domain for computation speed
  if( D==1 ){
    ##### Compute bootstrap means
    bootMeans <- R%*%multiplier / N

    ### Estimate the variance from the sample
      if( method == "regular" ){
        data.sigma <- sqrt(matrixStats::rowVars(R))
      }else if( method == "t" ){
        bootSecMoments <- R^2 %*% multiplier^2 / N
        # We put an abs here to make sure that no NaNs are produced due to machine precision error.
        data.sigma <- sqrt((N / (N-1)) * abs(bootSecMoments - bootMeans^2))
      }

    ### Compute bootstrap distribution of the maximum
    distVec <- sqrt(N)*apply( abs( bootMeans/data.sigma ), 2, max)

  }else if( D>1 ){
    ##### Compute bootstrap means
    bootMeans = array( matrix(R, prod(dimR[-(D+1)]), N)%*%multiplier, dim=c(dimR[-(D+1)], Mboots) ) / N

    ### Estimate the variance from the sample
    if( method == "regular" ){
      data.sigma <- sqrt( apply( R, 1:D, var  ) )
      distVec    <- apply( array( as.vector(abs(bootMeans)) / as.vector(data.sigma), dim=c(dimR[-(D+1)], Mboots) ), D+1, max)  * sqrt(N)
    }else if( method == "t" ){
      bootSecMoments <- array( matrix(R^2, prod(dimR[-(D+1)]), N)%*%multiplier^2, dim=c(dimR[-(D+1)], Mboots) ) / N
      # We put an abs here to make sure that no NaNs are produced due to machine precision error.
      data.sigma <- sqrt((N / (N-1)) * abs(bootSecMoments - bootMeans^2))
      distVec    <- apply( abs(bootMeans)  / data.sigma, D+1, max)  * sqrt(N)
    }
  }

  ### Return quantile and bootstrap distribution
  return( list( z = distVec, q = quantile(distVec, 1-alpha, type=8)) )
}
