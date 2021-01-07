################################################################################################
####                                                                                        ####
####                      Miscellaneous useful functions                                    ####
####                                                                                        ####
################################################################################################
## required packages:
##
## included functions:
##   - scaleField
##   - SNRresiduals
##   - SkewResiduals
##   - KurtResiduals
##
################################################################################################
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
