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

#' Computes the SNR residuals for a sample of 1d functional data
#'
#' @param Y Matrix containing the data. Last dimension must indicate different samples.
#' @param bias Boolean TRUE means variance estimated using 1/N, FALSE estimates it with 1/(N-1)
#' @return Scale field
#' @export
residualsSNR <- function(Y, bias=TRUE){
  # Get the dimension of the field
  dimY   = dim(Y);
  N      = dimY[length(dimY)];

  # Compute the factor correction for biased variance
  if(length(dimY)<=2){
      # Compute the sample mean and the sample variance
      mY = rowMeans(Y);
      sd = sqrt(matrixStats::rowVars(Y));

      # Compute the standard residuals, the asymptotic variance and the SNR
      R        = sqrt(N/(N-1))*(Y-mY);
      SNR      = mY/sd;
      asymptsd = sqrt( 1+SNR^2/2 );

      # Compute the modified residuals random variables with asymptotically the correct distribution
      res = ( R/sd - mY/(2*sd)*((R/sd)^2-1) );
  }else{
      ### Pointwise sample means and Pointwise sample variances
      mY   = array( rep(apply( Y, 1:D, mean ),N), dim = c(dimY[1:D], N) );
      sdY  = array( rep(sqrt(apply( Y, 1:D, var )),N), dim = c(dimY[1:D], N) );

      # Compute the standard residuals and the asymptotic variance
      R        = sqrt(N/(N-1))*(Y-mY);
      SNR      = mY / sdY;
      asymptsd = sqrt(1+SNR^2/2 );

      # Compute the modified residuals with asymptotically the correct distribution to estimate the LKCs
      res = ( R/sdY - mY/(2*sdY)*((R/sdY)^2-1) );
  }
  return( list( SNR=SNR, res=res, asymptsd=asymptsd) )

}
