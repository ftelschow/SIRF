################################################################################################
##                                                                                          ####
##               Functions to compute simultaneous confidence bands for functional data     ####
##                                                                                          ####
################################################################################################
## required packages:
##
## included functions:
##   - scb_mean (tested)
##   - scb_meandiff (tested)
##   - scb_glm (not included yet)
##   - scb_SNR (tested 1D)
##
################################################################################################
#' Computes simultaneous confidence bands for the mean of a sample from a one dimensional functional signal plus noise model. It is possible to choose between different estimators for the quantile.
#'
#' @param Y Array of dimension K_1 x ... x K_d x N containing N-realizations of a Gaussian random field over a d-dimensional domain.
#' @param level Numeric the targeted covering probability. Must be strictly between 0 and 1.
#' @param method Stringspecifying the method to construt the scb, i.e. estimatate the quantile. Current options are "tGKF", "GKF", "Bootstrap", "Bootstrapt", "MultiplierBootstrap", "MultiplierBootstrapt". Default value is "tGKF".
#' @param param_method list containing the parameters for 'method'. The list must contain the following elements, otherwise default values are set:
#'  \itemize{
#'   \item For method either "tGKF" or "GKF": L0 (integer) the Euler characteristic of the domain of the random fields Y [default is 1], LKC_estim (function) a function estimating the Lipschitz-Killing curvatures from the normalized residuals [default is LKC_estim_direct].
#'   \item For method "Bootstrap", "Bootstrapt", "MultiplierBootstrap", "MultiplierBootstrapt": Mboots (positiv integer) the amount of bootstrap replicates [default 5e3].
#' }
#'
#' @return list with elements
#'  \itemize{
#'   \item hatmean pointwise sample mean
#'   \item scb list containing the upper and lower bounds of the simultaneous confidence band
#'   \item level targeted covering probability
#'   \item q quantile of the maximum of the residual field
#' }
#' @export
scb_mean <- function( Y, level=.95, method="tGKF", param_method=NULL ){
  ###### Check user input
  ### Check input Y
  if(!is.array(Y)){
    stop("Y must be an array containing the realisations of the functional data, i.e. an array with last dimension enumerating the realisations.")
  }

  ### Check input level
  if( is.numeric(level) ){
    if( level>=1 | level <=0 ){
      stop("The input 'level' needs to be strictly between 0 and 1.")
    }
  }else{stop("The input 'level' needs to be strictly between 0 and 1.")}

  ### Check input method
  if( is.character(method) ){
    if( !(method%in%c("tGKF", "GKF", "Bootstrap", "Bootstrapt", "MultiplierBootstrap", "MultiplierBootstrapt")) ){
      stop("Choose a valid option from the available quantile approximations. Please, check the help page.")
    }
  }else{stop("Choose a valid option from the available quantile approximations. Please, check the help page.")}

  ### Check input param_method
  if( is.list(param_method) | is.null(param_method) ){
    ## Check input for methods "tGKF" and "GKF" and put the default values, if necessary.
    if( method%in%c("tGKF", "GKF") ){
          ## Check the Euler Characterisitc input L0
          if( is.null(param_method$L0) ){
            param_method$L0 = 1
          }else{
            if( is.numeric(param_method$L0) ){
              if( param_method$L0%%1 !=0 ){
                stop("The element param_method$L0 must be an integer.")
              }
            }else{
              stop("The element param_method$L0 must be an integer.")
            }
          }
          ## Check the LKC_estim input
          if( is.null(param_method$LKC_estim) ){
            param_method$LKC_estim = LKC_estim_direct
          }else{
            if( !is.function(param_method$LKC_estim) ){
              stop("The element param_method$LKC_estim must be a function computing the Lipschitz Killing curvatures from the normed residuals.")
            }
          }
    }
    ## Check input for methods "Bootstrap", "Bootstrapt", "MultiplierBootstrap", "MultiplierBootstrapt" and put the default values, if necessary.
    if( method%in%c("Bootstrap", "Bootstrapt", "MultiplierBootstrap", "MultiplierBootstrapt") ){
          if( is.null(param_method$Mboots) ){
            param_method$Mboots = 5e3
          }else{
            if( is.numeric(param_method$Mboots) ){
              if( param_method$Mboots%%1 !=0 | param_method$Mboots < 1  ){
                stop("The element param_method$Mboots must be a positiv integer.")
              }
            }else{
              stop("The element param_method$Mboots must be a positiv integer.")
            }
          }
    }
  }else{
    stop("'param_method' must be a list containing the elements specified in the help page or NULL, if you want to use the standard options of the method.")
  }

  ###### Compute useful constants
  dimY = dim(Y);
  ### dimension of domain
  D    = length(dimY)-1;
  ### get number of sample curves
  N    = dimY[length(dimY)]

  ###### Compute the necessary statistics from the data for residuals
  ### Pointwise sample means
  mY   = array( rep(apply( Y, 1:D, mean ),N), dim = c(dimY[1:D], N) );
  ### Pointwise sample variances
  sd2Y = array( rep(apply( Y, 1:D, var ),N), dim = c(dimY[1:D], N) );

  ###### Estimate the quantile of the maximum of the absolute value of the limiting Gaussian process of the CLT for the mean
  if( method=="tGKF" ){
    ### Compute the residuals to estimate the LKCs
    R   = (Y-mY) / sqrt( sd2Y );
    ### Estimate the LKCs
    LKC = c(param_method$L0, param_method$LKC_estim(R));
    q   = GKFquantileApprox( alpha = (1-level)/2, LKC, field="t", df=N-1 );
  }else if( method=="GKF" ){
    ### Compute the residuals to estimate the LKCs
    R = (Y-mY) / sqrt( sd2Y );
    ### Estimate the LKCs
    LKC = c(param_method$L0, param_method$LKC_estim(R));
    q   = GKFquantileApprox( alpha = (1-level)/2, LKC, field="Gauss", df=1 );
  }else if( method=="Bootstrap" ){
    ### Compute the residuals to estimate the LKCs
    R = Y-mY
    ### Estimate the quantile
    q <- NonParametricBootstrap( R, params=list( Mboots = param_method$Mboots, alpha = 1-level, method="regular" ) )$q
  }else if( method=="Bootstrapt" ){
    ### Compute the residuals to estimate the LKCs
    R = Y-mY
    ### Estimate the quantile
    q <- NonParametricBootstrap( R, params=list( Mboots = param_method$Mboots, alpha = 1-level, method="t" ) )$q
  }else if( method=="MultiplierBootstrap" ){
    ### Compute the residuals to estimate the LKCs
    R = sqrt(N / (N-1))*(Y-mY)
    ### Estimate the quantile
    q <- MultiplierBootstrap( R, params=list( Mboots = param_method$Mboots, alpha = 1-level, method="regular" ) )$q
  }else{
    ### Compute the residuals to estimate the LKCs
    R = sqrt(N / (N-1))*(Y-mY)
    ### Estimate the quantile
    q <- MultiplierBootstrap( R, params=list( Mboots = param_method$Mboots, alpha = 1-level, method="t") )$q
  }

  ###### compute the simultaneous confidence bands
  ### Pointwise sample means
  mY   = apply( Y, 1:D, mean );
  ### Pointwise sample variances
  sd2Y = apply( Y, 1:D, var );

  scb    <- list()
  scb$lo <- mY - q *  sqrt( sd2Y ) / sqrt( N )
  scb$up <- mY + q *  sqrt( sd2Y ) / sqrt( N )

  ### Return a list containing estimate of the mean, SCBs etc.
  list( hatmean = mY, scb = scb, level = level, q = q, res = R  )
}

#' Computes simultaneous confidence bands for the difference of the means of two samples from one dimensional functional signal plus noise models. It is possible to choose between different estimators for the quantile.
#'
#' @param Y1 Array dimension K_1 x ... x K_d x N containing N-realizations of a Gaussian random field over a d-dimensional domain. This is sample 1.
#' @param Y2 Array dimension K_1 x ... x K_d x N containing N-realizations of a Gaussian random field over a d-dimensional domain. This is sample 2.
#' @param level Numeric containing the level of confidence for the constructed SCBs.
#' @param method String the name of the method used for quantile estimation. Possibilities are "tGKF" and "GKF". Default value is "tGKF".
#' @return LKC Vector containing the Lipschitz killing curvatures of dimension greater than 1. Note that the 0-th LKC must be calculated seperately as the Euler characteristic of the domain.
#' @export
scb_meandiff <- function( Y1, Y2, level = .95, method="tGKF", param_method=NULL ){
  ###### Check user input
  ## Check Y1
  if(!is.array(Y1)){
    stop("Y1 must be an array containing the realisations of the functional data.")
  }

  ## Check Y2
  if(!is.array(Y2)){
    stop("Y2 must be an array containing the realisations of the functional data.")
  }

  ## Check level
  if( is.numeric(level) ){
    if( level>=1 | level <=0 ){
      stop("The input 'level' needs to be a number strictly between 0 and 1.")
    }
  }else{
    stop("The input 'level' needs to be a number strictly between 0 and 1.")
  }

  ## Check method
  if( is.character(method) ){
    if( !(method%in%c("tGKF", "GKF", "MultiplierBootstrap", "MultiplierBootstrapt")) ){
      stop("Choose a valid option from the available quantile approximations. Check help page.")
    }
  }else{
    stop("Choose a valid option from the available quantile approximations. Check help page.")
  }

  ## Check input param_method
  if( is.list(param_method) | is.null(param_method) ){
      ## Check input for methods "tGKF" and "GKF" and put the default values, if necessary.
      if( method%in%c("tGKF", "GKF") ){
          ## Check the Euler Characterisitc input L0
          if( is.null(param_method$L0) ){
              param_method$L0 = 1
          }else{
              if( is.numeric(param_method$L0) ){
                  if( param_method$L0%%1 !=0 ){
                      stop("The element param_method$L0 must be an integer.")
                  }
              }else{
                  stop("The element param_method$L0 must be an integer.")
              }
          }
          ## Check the LKC_estim input
          if( is.null(param_method$LKC_estim) ){
              param_method$LKC_estim = LKC_estim_direct
          }else{
              if( !is.function(param_method$LKC_estim) ){
              stop("The element param_method$LKC_estim must be a function computing the Lipschitz Killing curvatures from the normed residuals.")
              }
          }
      }else if( method%in%c("MultiplierBootstrap", "MultiplierBootstrapt") ){
            ## Check the Msim input
            if( is.null(param_method$Msim) ){
              param_method$Msim = 5e3
            }else{
              if( is.numeric(param_method$Msim) ){
                if( param_method$Msim%%1 !=0 | param_method$Msim<=0 ){
                  stop("The element param_method$Msim must be a positiv integer.")
                }
              }else{
                stop("The element param_method$Msim must be a positiv integer.")
              }
            }
      }else{
          stop("'param_method' must be a list containing the elements specified in the help page or NULL, if you want to use the standard options of the method.")
      }
  }

  ###### Get constants from the input
  ## Dimensions of the data
  dimY1 = dim(Y1);
  dimY2 = dim(Y2);
  ## Dimension of domain of the field
  D = length(dimY1)-1;
  ## Number of sample curves
  N1 = dimY1[D+1];
  N2 = dimY2[D+1];
  c = N2/N1;

  ###### Compute the necessary statistics from the data
  ### Pointwise sample means
  mY1   = array( rep(apply( Y1, 1:D, mean ), N1), dim = c(dimY1[1:D], N1) );
  mY2   = array( rep(apply( Y2, 1:D, mean ), N2), dim = c(dimY2[1:D], N2) );
  ### Pointwise sample variances
  sd2Y1 = apply( Y1, 1:D, var );
  sd2Y2 = apply( Y2, 1:D, var );

  ###### Compute the residuals to estimate the LKCs
  R1 = (Y1 - mY1) / sqrt( (1+c)*array( rep(sd2Y1,N1), dim = c(dimY1[1:D], N1) ) + (1+c^{-1})*array( rep(sd2Y2,N1), dim = c(dimY1[1:D], N1) ) ) * sqrt(1+c);
  R2 = (Y2 - mY2) / sqrt( (1+c)*array( rep(sd2Y1,N2), dim = c(dimY1[1:D], N2) ) + (1+c^{-1})*array( rep(sd2Y2,N2), dim = c(dimY1[1:D], N2) ) ) * sqrt(1+c^{-1});

  ###### Estimate the quantile of the maximum of the limiting Gaussian process
  if(method=="tGKF"){
    ### Estimate the LKCs
    LKC = c(param_method$L0, param_method$LKC_estim( R=abind::abind(R1, R2, along=D+1), subdiv=c(N1,N1+N2) ));
    q   = GKFquantileApprox( alpha = (1-level)/2, LKC, field="t", df=N1+N2-2 );
  }else if(method=="GKF"){
    ### Estimate the LKCs
    LKC = c(param_method$L0, param_method$LKC_estim( abind::abind(R1, R2, along=D+1), c(N1,N1+N2) ));
    q   = GKFquantileApprox( alpha = (1-level)/2, LKC, field="Gauss", df=1 );
  }else if( method=="MultiplierBootstrap" ){
    ### Modify the residuals
    R1 = sqrt(N1 / (N1-1))*R1
    R2 = sqrt(N2 / (N2-1))*R2
    ### Estimate the quantile
    q <- MultiplierBootstrap( R1-R2, params=list( Mboots = param_method$Mboots, alpha = 1-level, method="regular" ) )$q
  }else{
    ### Modify the residuals
    R1 = sqrt(N1 / (N1-1))*R1
    R2 = sqrt(N2 / (N2-1))*R2
    ### Estimate the quantile
    q <- MultiplierBootstrap( R1-R2, params=list( Mboots = param_method$Mboots, alpha = 1-level, method="t" ) )$q
  }

  ###### compute the simultaneous confidence bands
  ## Pointwise sample means
  mY1 = apply( Y1, 1:D, mean );
  mY2 = apply( Y2, 1:D, mean );
  ## Pointwise sample variances
  sd2Y1 = apply( Y1, 1:D, var );
  sd2Y2 = apply( Y2, 1:D, var );
  ## SCB computation
  scb    <- list()
  scb$lo <- mY1-mY2 - q *  sqrt( (1+c)*sd2Y1 + (1+c^{-1})*sd2Y2 ) / sqrt( N1+N2 )#/ sqrt( N1+N2-2 )
  scb$up <- mY1-mY2 + q *  sqrt( (1+c)*sd2Y1 + (1+c^{-1})*sd2Y2 ) / sqrt( N1+N2 )#/ sqrt( N1+N2-2 )

  ### Return a list containing estimate of the mean, SCBs etc.
  list( hatmean = mY1-mY2, scb = scb, level = level, q = q, res1 = R1, res2 = R2  )
}


#' Computes simultaneous confidence bands for a contrast in a functional linear model. It is possible to choose between different estimators for the quantile.
#'
#' @param Y Array (dimension K_1 x ... x K_d x N containing N-realizations of a Gaussian random field over a d-dimensional domain.)
#' @param X Matrix (dimensions N x P. It is the design matrix of the pointwise functional linear model.)
#' @param c Vector (of length P defining the contrast of interest.)
#' @param level Numeric (the targeted covering probability. Must be strictly between 0 and 1.)
#' @param method String (the name of the method used for quantile estimation. Possibilities are "tGKF", "GKF", "Bootstrap", "Bootstrapt", "MultiplierBootstrap", "MultiplierBootstrapt". Default value is "tGKF".)
#' @param Mboots Numeric (the number of bootstrap replicats used, if the quantile is estimated using a bootstrap method. Default value is 5000.)
#' @return list with elements
#'  \itemize{
#'   \item hatmean pointwise sample mean
#'   \item scb list containing the upper and lower bounds of the simultaneous confidence band
#'   \item level targeted covering probability
#'   \item q quantile of the maximum of the residual field
#' }
#' @export
SCBglm <- function( x=seq(0,1,length.out=100), y, X=NULL, c=c(1,-1), xlim=c(0,1), level=.95, quantile.method="tGKF", Mboots=1000, evalN=ifelse( is.vector(x), 4*length(x), 4*length(x[1,]) ), LinWeights=NULL, bw=NULL, kernel=locpol::gaussK, sigma.boots="boots", quantile=NULL ){

  ### evaluation grid
  xeval = seq( xlim[1], xlim[2], length.out = evalN )
  ### get number of sample curves
  nSamp <- dim(y)[2]

  ###### Smooth data, if neccessary
  if( !is.null(LinWeights) ){
    ysmooth <- LinWeights %*% y
  }else if( is.null(LinWeights) & !is.null(bw) ){
    LinWeights  <- locpol::locLinWeightsC( x = x, xeval = xeval, bw = bw, kernel = kernel)$locWeig
    ysmooth     <- LinWeights %*% y
  }else{
    ysmooth     <- y
    xeval       <- x
  }
  ###### END Smooth data

  ### BLUS estimator and variance of the design
  XTX       = solve( t(X)%*%(X))
  designVar = as.vector(sqrt( t(c)%*%XTX%*%c ))

  ### Compute residuals from linear model
  hatbeta   = ysmooth%*%t(XTX%*%t(X))
  Residuals = ysmooth-hatbeta%*%t(X)
  hatsigma  = sqrt(matrixStats::rowVars(Residuals) * (nSamp-1)/nSamp)

  ### quantile estimation of maximum of t-field using GKF for t-fields
  q.alpha      <- gaussianKinematicFormulaT( A = Residuals, nu = nSamp-length(c), alpha = 1-level, center=FALSE, normalize=TRUE )$c.alpha

  ### construct confidence bands
  chatbeta = as.vector(hatbeta%*%c)

  scb  <- chatbeta + cbind( -q.alpha *  hatsigma*designVar * sqrt( nSamp ) / sqrt( nSamp-length(hatbeta) ), q.alpha *  hatsigma*designVar * sqrt( nSamp ) / sqrt( nSamp-length(hatbeta) ) )

  if( is.null(LinWeights) ){
    ### interpolate confidence bands and mean to get the values on output grid
    xeval   <- seq( xlim[1], xlim[2], length.out = evalN )
    chatbeta <- approx( x, chatbeta, xout=xeval )$y
    scb  <- apply( scb, 2, function(y) approx( x, y, xout=xeval )$y )

    retList <- list( x=x, y=y, hatmean=chatbeta, level=level, q=q.alpha, xeval=xeval, scb=scb, res=R )
  }else{
    retList <- list( x=x, y=y, LinWeights=LinWeights, hatmean=chatbeta, level=level, q=q.alpha, xeval=xeval, scb=scb, res=Residuals )
  }
  ### Return a list containing estimate of the mean, SCBs etc.
  retList
}


#' Computes simultaneous confidence bands for the SNR of a sample from a one dimensional functional signal plus noise model. It is possible to choose between different estimators for the quantile.
#'
#' @param Y Array of dimension K_1 x ... x K_d x N containing N-realizations of a Gaussian random field over a d-dimensional domain.
#' @param level Numeric the targeted covering probability. Must be strictly between 0 and 1.
#' @param method Stringspecifying the method to construt the scb, i.e. estimatate the quantile. Current options are "tGKF", "GKF", "Bootstrap", "Bootstrapt", "MultiplierBootstrap", "MultiplierBootstrapt". Default value is "tGKF".
#' @param param_method list containing the parameters for 'method'. The list must contain the following elements, otherwise default values are set:
#'  \itemize{
#'   \item For method either "tGKF" or "GKF": L0 (integer) the Euler characteristic of the domain of the random fields Y [default is 1], LKC_estim (function) a function estimating the Lipschitz-Killing curvatures from the normalized residuals [default is LKC_estim_direct].
#'   \item For method "Bootstrap", "Bootstrapt", "MultiplierBootstrap", "MultiplierBootstrapt": Mboots (positiv integer) the amount of bootstrap replicates [default 5e3].
#' }
#'
#' @return list with elements
#'  \itemize{
#'   \item hatmean pointwise sample mean
#'   \item scb list containing the upper and lower bounds of the simultaneous confidence band
#'   \item level targeted covering probability
#'   \item q quantile of the maximum of the residual field
#'   \item res residual field
#' }
#' @export
scb_SNR <- function( Y, level=.95, method="GKF", param_method=NULL ){
  ###### Check user input
  ### Check input Y
  if(!is.array(Y)){
    stop("Y must be an array containing the realisations of the functional data, i.e. an array with last dimension enumerating the realisations.")
  }

  ### Check input level
  if( is.numeric(level) ){
    if( level>=1 | level <=0 ){
      stop("The input 'level' needs to be strictly between 0 and 1.")
    }
  }else{stop("The input 'level' needs to be strictly between 0 and 1.")}

  ### Check input method
  if( is.character(method) ){
    if( !(method%in%c("tGKF", "GKF", "Bootstrap", "Bootstrapt", "MultiplierBootstrap", "MultiplierBootstrapt")) ){
      stop("Choose a valid option from the available quantile approximations. Please, check the help page.")
    }
  }else{stop("Choose a valid option from the available quantile approximations. Please, check the help page.")}

  ### Check input param_method
  if( is.list(param_method) | is.null(param_method) ){
    ## Check input for methods "tGKF" and "GKF" and put the default values, if necessary.
    if( method%in%c("tGKF", "GKF") ){
      ## Check the Euler Characterisitc input L0
      if( is.null(param_method$L0) ){
        param_method$L0 = 1
      }else{
        if( is.numeric(param_method$L0) ){
          if( param_method$L0%%1 !=0 ){
            stop("The element param_method$L0 must be an integer.")
          }
        }else{
          stop("The element param_method$L0 must be an integer.")
        }
      }
      ## Check the LKC_estim input
      if( is.null(param_method$LKC_estim) ){
        param_method$LKC_estim = LKC_estim_direct
      }else{
        if( !is.function(param_method$LKC_estim) ){
          stop("The element param_method$LKC_estim must be a function computing the Lipschitz Killing curvatures from the normed residuals.")
        }
      }
    }
    ## Check input for methods "Bootstrap", "Bootstrapt", "MultiplierBootstrap", "MultiplierBootstrapt" and put the default values, if necessary.
    if( method%in%c("Bootstrap", "Bootstrapt", "MultiplierBootstrap", "MultiplierBootstrapt") ){
      if( is.null(param_method$Mboots) ){
        param_method$Mboots = 5e3
      }else{
        if( is.numeric(param_method$Mboots) ){
          if( param_method$Mboots%%1 !=0 | param_method$Mboots < 1  ){
            stop("The element param_method$Mboots must be a positiv integer.")
          }
        }else{
          stop("The element param_method$Mboots must be a positiv integer.")
        }
      }
    }
  }else{
    stop("'param_method' must be a list containing the elements specified in the help page or NULL, if you want to use the standard options of the method.")
  }

  ###### Compute useful constants
  dimY = dim(Y);
  ### dimension of domain
  D    = length(dimY)-1;
  ### get number of sample curves
  N    = dimY[length(dimY)]

  ###### Compute the necessary statistics from the data for residuals
  ### Pointwise sample means
  mY   = array( rep(apply( Y, 1:D, mean ),N), dim = c(dimY[1:D], N) );
  ### Pointwise sample variances
  sdY  = array( rep(sqrt(apply( Y, 1:D, var )),N), dim = c(dimY[1:D], N) );

  ###### Estimate the quantile of the maximum of the absolute value of the limiting Gaussian process of the CLT for the mean
  if( method=="tGKF" ){
    # Compute the modified residuals with asymptotically the correct distribution to estimate the LKCs
    R = ( (Y-mY)/sdY - mY/(2*sdY)*(((Y-mY)/sdY)^2-1) ) / sqrt(1+mY^2/(2*sdY^2) );
    ### Estimate the LKCs
    LKC = c(param_method$L0, param_method$LKC_estim(R));
    q   = GKFquantileApprox( alpha = (1-level)/2, LKC, field="t", df=N-1 );
  }else if( method=="GKF" ){
    # Compute the modified residuals with asymptotically the correct distribution to estimate the LKCs
    R = ( (Y-mY)/sdY - mY/(2*sdY)*(((Y-mY)/sdY)^2-1) ) / sqrt(1+mY^2/(2*sdY^2) );
    ### Estimate the LKCs
    LKC = c(param_method$L0, param_method$LKC_estim(R));
    q   = GKFquantileApprox( alpha = (1-level)/2, LKC, field="Gauss", df=1 );
  }else if( method=="Bootstrap" ){
    # Compute the modified residuals with asymptotically the correct distribution to estimate the LKCs
    R = ( (Y-mY)/sdY - mY/(2*sdY)*(((Y-mY)/sdY)^2-1) );
    ### Estimate the quantile
    q <- NonParametricBootstrap( R, params=list( Mboots = param_method$Mboots, alpha = 1-level, method="regular" ) )$q
  }else if( method=="Bootstrapt" ){
    ### Compute the residuals to estimate the LKCs
    # Compute the modified residuals with asymptotically the correct distribution to estimate the LKCs
    R = ( (Y-mY)/sdY - mY/(2*sdY)*(((Y-mY)/sdY)^2-1) );
    ### Estimate the quantile
    q <- NonParametricBootstrap( R, params=list( Mboots = param_method$Mboots, alpha = 1-level, method="t" ) )$q
  }else if( method=="MultiplierBootstrap" ){
    # Compute the modified residuals with asymptotically the correct distribution to estimate the LKCs
    R = ( (Y-mY)/sdY - mY/(2*sdY)*(((Y-mY)/sdY)^2-1) ) / sqrt(1+mY^2/(2*sdY^2) );
    ### Estimate the quantile
    q <- MultiplierBootstrap( R, params=list( Mboots = param_method$Mboots, alpha = 1-level, method="regular" ) )$q
  }else{
    # Compute the modified residuals with asymptotically the correct distribution to estimate the LKCs
    R = ( (Y-mY)/sdY - mY/(2*sdY)*(((Y-mY)/sdY)^2-1) ) / sqrt(1+mY^2/(2*sdY^2) );
    ### Estimate the quantile
    q <- MultiplierBootstrap( R, params=list( Mboots = param_method$Mboots, alpha = 1-level, method="t") )$q
  }

  ###### compute the simultaneous confidence bands
  ### Pointwise sample means
  mY   = apply( Y, 1:D, mean );
  ### Pointwise sample variances
  sdY = sqrt(apply( Y, 1:D, var ));
  ### Pointwise estimate of the SNR and its asymptotic variance
  nu        = mY / sdY;
  sd_nuasym = sqrt(1+mY^2/(2*sdY^2) );

  scb    <- list()
  scb$lo <- nu - q *  sd_nuasym / sqrt( N )
  scb$up <- nu + q *  sd_nuasym / sqrt( N )

  ### Return a list containing estimate of the mean, SCBs etc.
  list( hatnu = nu, scb = scb, level = level, q = q, res=R )
}
