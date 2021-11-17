#------------------------------------------------------------------------------#
#                                                                              #
#     Miscalaneous functions                                                   #
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#  - scaleField
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline NOT included
#------------------------------------------------------------------------------#
#' This functions computes the standard
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
se_Gaussian <- function(transformation, N,...){
  if(transformation == "cohensd"){# sqrt(N)\hat d ~ non-central t with nc = sqrt(N)d, df = N-1
    df = N - 1
    nc = sqrt(N) * hatd

    if( df < 320){
      fac <- sqrt(df / 2) * gamma((df - 1) / 2) / gamma(df / 2)
    }else{
      fac <- 1 / (1 - 3 / (4 * df -1 ))
    }

    se <- sqrt(( df * (1 + nc^2) / (df - 2) - nc^2 * fac^2 )) / sqrt(N)
  }
  else if(transformation == "skewness"){
    se <- sqrt(6 * (N - 2) / (N + 1) / (N + 3))
  }else if(transformation == "kurtosis"){
    se <- sqrt(24 * N * (N - 2) * (N - 3) / (N + 1)^2 / (N + 3) / (N + 5))
  }
  return(se)
}

#' This functions computes the bias of certain transformations under Gaussianity
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
bias_Gaussian <- function(transformation, N,...){
  if(transformation == "cohensd"){# sqrt(N)\hat d ~ non-central t with nc = sqrt(N)d, df = N-1
    nu = N-1

    if( df < 320){
      fac <- sqrt(df / 2) * gamma((df - 1) / 2) / gamma(df / 2)
    }else{
      fac <- 1 / (1 - 3 / (4 * df -1 ))
    }

    bias <- hatd * fac
  }
  else if(transformation == "skewness"){
    bias <- 0
  }else if(transformation == "kurtosis"){
    bias <- -6 / (N + 1)
  }
  return(bias)
}

#' Computes the scale field for a sample of 1d functional data
#'
#' @param Y Matrix containing the data. Last dimension must indicate different samples.
#' @param xeval Vector grid for evaluation of the smoothed field
#' @param h Vector grid for the bandwidths for local linear smoothing or spacing between the knots in B-spline smoothing
#' @param method String Either 'loclin' for local linear smoothing or 'smoothingspline' for Bspline smoothing
#' @param kernel Function the kernel used for smoothing the data
#' @param Weights Matrix containing a smoothing matrix. If specified all other parameters for smoothing are ignored.
#' @return Scale field
#' @export
scaleField <- function( Y, x=seq(0,1,length.out=nrow(Y)), xeval=seq(0,1,length.out=2*nrow(Y)), h = seq(x[2],0.1,length.out=20), method="loclin", kernel=locpol::gaussK, Weights=NULL ){
  ########### Check User Input
  ### Check Y
  if( !(is.matrix(Y) | is.vector(Y)) ){
    stop("Y must be either a vector or a matrix.")
  }

  ##### If Y is a vector change it to a matrix
  if( is.vector(Y) ){
    Y = matrix(Y, length(Y),1)
  }

  ##### Divide cases by dimension
  if( length(dim(Y))==2 ){
    #### Smooth the data with your chosen favourite method (default is local linear)
    if( !is.null(Weights) ){
      H = dim(Weights)[3]
      #### Initialize the scale field as an array
      ScaleField    = array(NA, c( dim(Weights)[1], H, dim(Y)[2] ) )
      #### Compute smoothed data
      for(j in 1:H){
        ScaleField[,j,] <- Weights[,,j] %*% Y
      }
    }else if( method=="loclin" ){
      #### Initialize the scale field as an array
      ScaleField    = array(NA, c( length(xeval), length(h), dim(Y)[2] ) )

      SmoothWeights = matrix(NA, length(xeval), length(x))
      for( j in 1:length(h) ){
        #### Compute weights for local linear regression
        SmoothWeights   <- as.matrix( locpol::locLinWeightsC( x = x, xeval = xeval, bw = h[j], kernel = kernel )$locWeig )

        #### Compute smoothed data
        ScaleField[,j,] <- SmoothWeights %*% Y
      }
    }else if( method=="smoothingspline" ){
      #### Initialize the scale field as an array
      ScaleField    = array(NA, c( length(xeval), length(h), dim(Y)[2] ) )
      for( j in 1:length(h) ){
        #### Compute smoothed data
        ScaleField[,j,] <- vapply( 1:dim(Y)[2], function(n) smooth.spline(x, y = Y[,n], lambda=h[j], all.knots = TRUE  )$y, FUN.VALUE=rep(0,length(xeval)) )
      }
    }else{
      stop( "Please, specify as method either 'loclin' or 'bsplines'" )
    }
  }else{
    stop( "Not implemented for data of dimension 2 or greater" )
  }
  return( ScaleField )
}
