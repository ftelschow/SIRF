################################################################################################
##                                                                                            ##
##  This file contains EC densities and Gaussian kinematic formulas, as well as quantile      ##
##  estimators for the maximum of random fields using the Euler characteristic heuristic      ##
##                                                                                            ##
################################################################################################
#
# Included functions:
#   - EC_density
#   - GKFquantileApprox
#   - scb_glm (not included yet)
#   - scb_SNR (not included yet)
#
################################################################################################
#' Computes the d-th EC density of different fields. Currently only "t" is implemented.
#'
#' @param x Numeric location x at which the EC density is to be evaluated
#' @param d Numeric d-th EX density
#' @param field String type of field for EC densities, currently only "t" (Student-t) and "Gauss" (Gaussian) are supported
#' @param df Numeric degree of freedom of t-field
#' @return value of 2d-EC density of a t-field of degree df at location x.
#' @export
ECdensity <- function( x, d, field="t", df=1 ){
  if( field=="t" ){
      if(d==0){
         1 - pt( x, df = df )
      }else if(d==1){
        (2*pi)^(-1) * ( 1 + x^2/df )^( -(df-1)/2 )
      }else if(d==2){
        (2*pi)^(-3/2) * gamma( (df+1)/2 ) / gamma( df/2 ) / sqrt(df) * sqrt(2) * x * ( 1 + x^2/df )^( -(df-1)/2 )
      }else if(d==3){
        (2*pi)^(-2) * ( (df-1)/df*x^2 - 1 ) * ( 1 + x^2/df )^( -(df-1)/2 )
      }else{
        stop("Error: d must be smaller then 3. Higher dimensions are not yet implemented.")
      }
  }else if( field=="Gauss" ){
      constFunc = (2*pi)^(-(d+1)/2)*exp(-x^2/2)
      if(d==0){
        1 - pnorm( x )
      }else if(d==1){
        constFunc * 1
      }else if(d==2){
        constFunc * (2*x)
      }else if(d==3){
        constFunc * (4*x^2-2)
      }else{
        stop("Error: d must be smaller then 4. Higher dimensions are not yet implemented.")
      }
  }else{
    stop("Error: Input a valid field type.")
  }
}

#' Approximates the upper tail probabilities of different fields using the corresponding Gaussian
#' kinematic formulas.
#' It approximates the tail probability of the absolute value by twice the
#' usual tail probability of the process.
#'
#' @param alpha Numeric upper tail probabiltiy, which needs to be approximated.
#' @param LKC Vector estimated LKCs of the UGRFs of the field.
#' @param field String type of field.
#' @param df degrees of freedom.
#' @return Numeric the approximated 1-alpha quantile.
#' @export
GKFquantileApprox <- function( alpha = 0.05, LKC, field="t", df ){
  ###### Approximate the tail distribution using the GKF and EECH
  if( length(LKC) == 2 ){
    ### 1D case
    tailProb <- function(u){
            LKC[1] * ECdensity( x=u, d=0, field=field, df=df ) + LKC[2] * ECdensity( x=u, d=1, field=field, df=df ) - alpha
    }
  }else if( length(LKC) == 3 ){
    ### 2D case
    tailProb <- function(u){
        LKC[1] * ECdensity( x=u, d=0, field=field, df=df ) + LKC[2] * ECdensity( x=u, d=1, field=field, df=df ) + LKC[3] * ECdensity( x=u, d=2, field=field, df=df ) - alpha
    }
  }else if( length(LKC) == 4 ){
    ### 2D case
    tailProb <- function(u){
        LKC[1] * ECdensity( x=u, d=0, field=field, df=df ) + LKC[2] * ECdensity( x=u, d=1, field=field, df=df ) + LKC[3] * ECdensity( x=u, d=2, field=field, df=df ) + LKC[4] * ECdensity( x=u, d=3, field=field, df=df )
      - alpha
    }
  }else{
    stop( "Error: Not yet implemented for data with a domain of dimension greater than 3!")
  }
  uniroot( tailProb, interval = c( 0, 50 ) )$root
}
