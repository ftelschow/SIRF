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
  }else{
    stat = params$stat
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
      N     = dim(A)
      meanY = rowMeans(A); sdY = sqrt(matrixStats::rowVars(A))

      boots.SNR = matrix(NA, N[1], Mboots)

      for(m in 1:Mboots){
        Yboots = Y[,sample(1:N[2], replace=TRUE)]

        boots.mean    <- rowMeans(Yboots) #[,m] <- rowMeans(Yboots)
        boots.sd      <- sqrt(matrixStats::rowVars(Yboots))#[,m]  <- sqrt(matrixStats::rowVars(Yboots))
        ### Estimate the variance from the sample if it is unknown
        if( method == "regular" ){
          boots.SNR[,m] <- sqrt(N[2])*( boots.mean/boots.sd - meanY/sdY ) / sqrt(1+meanY^2/2/sdY^2)
        }else if( method == "t" ){
          boots.SNR[,m] <- sqrt(N[2])*( boots.mean/boots.sd - meanY/sdY ) / sqrt(1+boots.mean^2/2/boots.sd^2)
        }
      }

      ### Compute the bootstrap realisation of the maximum of the T-field approximation
      distVec <- apply(abs(boots.SNR), 2, max )

    }else if( D>1 ){

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
MultiplierBootstrap <- function( R, Q=NULL, params=NULL ){
  ###### Check user input and put default values
  #### Check R
  if( !is.array(R) ){
    stop("'R' must be an array.")
  }
  #### Check Q
  if( !(is.array(Q) | is.null(Q)) ){
    stop("'Q' must be an array.")
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
  ## Check params$weights
  if( is.null(params$weights) ){
    weights = "rademacher"
  }else{
    if( is.character(params$weights) ){
      if( !(params$weights %in% c("gauss", "rademacher", "mammen")) ){
        stop("Please, specify a valid choice for 'params$weights'. Options are 'gauss' and 'rademacher'.")
      }else{
        weights = params$weights
      }
    }else{
      stop("Please, specify a valid choice for 'params$weights'. Options are 'gauss' and 'rademacher'.")
    }
  }

  ###### One sample case
  if(is.null(Q)){
      ###### Precompute useful constants
      # dimension of input
      dimR = dim(R)
      # number of samples
      N    = dimR[length(dimR)]
      D    = length(dimR)-1

      ##### Simulate multiplier weights
      if(weights=="gauss"){
        multiplier <- matrix(rnorm(N * Mboots), N, Mboots)
      }else if(weights=="rademacher"){
        multiplier <- matrix(sample(c(-1,1), N*Mboots, replace=T), N, Mboots)
      }else{
        multiplier <- matrix(sqrt(5)*rbinom( N*Mboots, 1, (sqrt(5)-1)/2/sqrt(5)) + (1-sqrt(5))/2, N, Mboots)
      }
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
  }else{
    ###### Precompute useful constants
    # dimension of input
    dimR = dim(R)
    dimQ = dim(Q)
    # number of samples
    N    = dimR[length(dimR)]
    M    = dimQ[length(dimQ)]
    c    = N/M
    D    = length(dimR)-1

    ##### Simulate multiplier weights
    if(weights=="gauss"){
      multiplierR <- matrix(rnorm(N * Mboots), N, Mboots)
      multiplierQ <- matrix(rnorm(M * Mboots), M, Mboots)
    }else{
      multiplierR <- matrix(sample(1:2, N*Mboots, replace=T)*2-3, N, Mboots)
      multiplierQ <- matrix(sample(1:2, M*Mboots, replace=T)*2-3, M, Mboots)
    }

    ##### Split the cases into data with 1D domain and higher dimensional domain for computation speed
    if( D==1 ){
      ##### Compute bootstrap means
      bootMeansR <- R%*%multiplierR / N
      bootMeansQ <- Q%*%multiplierQ / M

      ### Estimate the variance from the samples
      if( method == "regular" ){
        data.sigmaR <- matrixStats::rowVars(R)
        data.sigmaQ <- matrixStats::rowVars(Q)
      }else if( method == "t" ){
        bootSecMoments <- R^2 %*% multiplierR^2 / N
        data.varR      <- (N / (N-1)) * abs(bootSecMoments - bootMeansR^2)

        bootSecMoments <- Q^2 %*% multiplierQ^2 / M
        data.varQ     <- (M / (M-1)) * abs(bootSecMoments - bootMeansQ^2)
      }
      data.sigma  <- sqrt( (1+1/c)*data.varR + (1+c)*data.varQ )

      ### Compute bootstrap distribution of the maximum
      distVec <- sqrt(N+M)*apply( abs( (bootMeansR + bootMeansQ) /data.sigma ), 2, max)

    }else if( D>1 ){
      ##### Compute bootstrap means
      bootMeansR = array( matrix(R, prod(dimR[-(D+1)]), N)%*%multiplierR, dim=c(dimR[-(D+1)], Mboots) ) / N
      bootMeansQ = array( matrix(Q, prod(dimR[-(D+1)]), M)%*%multiplierQ, dim=c(dimQ[-(D+1)], Mboots) ) / M

      ### Estimate the variance from the sample
      if( method == "regular" ){
        data.varR <- apply( R, 1:D, var  )
        data.varQ <- apply( Q, 1:D, var  )
        data.sigma  <- sqrt( (1+1/c)*data.varR + (1+c)*data.varQ )

        distVec    <- apply( array( as.vector(abs(bootMeans)) / as.vector(data.sigma), dim=c(dimR[-(D+1)], Mboots) ), D+1, max)  * sqrt(N+M-2)
      }else if( method == "t" ){
        #### bootstrapped variance for R
        bootSecMoments <- array( matrix(R^2, prod(dimR[-(D+1)]), N)%*%multiplierR^2, dim=c(dimR[-(D+1)], Mboots) ) / N
        data.varR <- (N / (N-1)) * abs(bootSecMoments - bootMeansR^2)
        #### bootstrapped variance for Q
        bootSecMoments <- array( matrix(Q^2, prod(dimQ[-(D+1)]), N)%*%multiplierQ^2, dim=c(dimQ[-(D+1)], Mboots) ) / M
        data.varQ <- (M / (M-1)) * abs(bootSecMoments - bootMeansQ^2)
        #### bootstrapped variance of the asymptotic process
        data.sigma  <- sqrt( (1+1/c)*data.varR + (1+c)*data.varQ )

        ##### max statistic
        distVec    <- apply( abs(bootMeansR + bootMeansQ)  / data.sigma, D+1, max)  * sqrt(N+M)
      }
    }
  }

  ### Return quantile and bootstrap distribution
  return( list( z = distVec, q = quantile(distVec, 1-alpha, type=8)) )
}
