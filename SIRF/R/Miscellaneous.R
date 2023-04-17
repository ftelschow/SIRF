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
#' This functions computes the Wilson & Hilferty transform
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
WF_transform <- function(x, df = 2){
  m  <- 1 - 2 / 9 / df
  sd <- sqrt(2 / 9 / df)

  ((x / df)^(1/3) - m) / sd
}

#' This functions computes the Wilson & Hilferty transform
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
Z1_transform <- function(g1, N){
  if(is.infinite(N)){
    asinh(g1 * 1.216675) * 0.3355443
    asinh(g1 * 1.231789) / 3
  }else{
    m21     <- var_Gaussian("skewness", N)
    gamma21 <- gamma2_Gaussian("skewness", N)

    W2 <- sqrt(2 * gamma21 + 4) - 1
    W2 <- sqrt(2 * gamma21 - 2 ) - 1
    d  <- 1 / sqrt(log(sqrt(W2))) / sqrt(N)
    a  <- sqrt(2 / (W2 - 1))

    d * asinh(g1 / a / sqrt(m21))
  }
}

#' This functions computes the Wilson & Hilferty transform
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
Z1_invtransform <- function(x, N){
  m21     <- var_Gaussian("skewness", N)
  gamma21 <- gamma2_Gaussian("skewness", N)

  W2    <- sqrt(2 * gamma21 + 4) - 1
  d <- 1 / sqrt(log(sqrt(W2))) / sqrt(N)
  a <- sqrt(2 / (W2 - 1))

  sinh(x / d) * a * sqrt(m21)
}

#' This functions computes the Wilson & Hilferty transform
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
Z2_transform <- function(g2, N){
  sq3 <- function(x){
    ifelse(x > 0, x^(1/3), -(-x)^(1/3))
  }

  if(N == Inf){
    return(0.8164966 * (1 - sq3(1 / (1 + (g2 - 3) * 0.75)  ) ))
  }else{
    m12 <- mean_Gaussian("kurtosis", N)
    sd22 <- sqrt(var_Gaussian("kurtosis", N))
    gamma12 <- gamma1_Gaussian("kurtosis", N)

    A <- 6 + 8 / gamma12 * (2 / gamma12 + sqrt(1 + 4 / gamma12^2))
    d2 <- sqrt(9 * A / 2 / N)

    return(d2 * (1 - 2 / 9 / A - sq3((1 - 2 / A) / (1 + (g2 - m12) / sd22 * sqrt(2 / (A - 4))  ) )))
  }
}

#' This functions computes the Wilson & Hilferty transform
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
Z2_invtransform <- function(x, N){
  m12 <- mean_Gaussian("kurtosis", N)
  sd22 <- sqrt(var_Gaussian("kurtosis", N))
  gamma12 <- gamma1_Gaussian("kurtosis", N)

  A <- 6 + 8 / gamma12 * (2 / gamma12 + sqrt(1 + 4 / gamma12^2))
  d2 <- sqrt(9 * A / 2 / N)

  ((1 - 2 / A) / (-x / d2 + (1 - 2 / 9 / A))^3 - 1) * sd22 / sqrt(2 / (A - 4)) + m12
}

#' This functions computes the values of Ksquare transformation under Gaussianity
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
Ksquare_Gaussian <- function(N){
  #---- Skewness transformation
  m21     <- 6 * (N - 2) / (N + 1) / (N + 3)
  gamma21 <- 36 * (N - 7) * (N^2 + 2*N - 5) / ((N - 2) * (N + 5) * (N + 7) * (N + 9))

  W2    <- sqrt(2 * gamma21 + 4) - 1
  delta1 <- 1 / sqrt(log(sqrt(W2))) / sqrt(N)
  alpha <- sqrt(2 / (W2 - 1))

  #---- kurtosis transformation
  m12 <- 3 * (N - 1) / (N + 1)
  m22 <- 24 * N * (N - 2) * (N - 3) / (N + 1)^2 / (N + 3) / (N + 5)
  gamma12 <- 6 * (N^2 - 5*N + 2) / (N + 7) / (N + 9) * sqrt(6 * (N + 3) * (N + 5) / N / (N - 2) / (N - 3))

  A <- 6 + 8 / gamma12 * (2 / gamma12 + sqrt(1 + 4 / gamma12^2))
  delta2 <- sqrt(9 * A / 2 / N)


  #---- K2 transformation
  K2 <- function(x,y){
    delta2^2 * ((1 - 2 / 9 / A) - ( (1 - 2 / A) /
                                      (1 + (y - m12)  / sqrt(m22) * sqrt(2 / (A - 4))) )^(1/3) )^2 +
      delta1^2 * log( x /
                        sqrt(m21) / alpha +
                        sqrt(  x / m21 / alpha^2 + 1) )^2
  }

  K2(0,3)
}

#' This functions computes the true standard error depending on N for different
#' functionals of Gaussians
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
mean_Gaussian <- function(transformation, N, hatd = NULL){
  if(transformation == "cohensd"){
    m <- NaN
  }
  else if(transformation == "skewness"){
    m <- 0
  }else if(transformation == "kurtosis"){
    m <- 3 * (N - 1) / (N + 1)
  }
  return(m)
}

#' This functions computes the true standard error depending on N for different
#' functionals of Gaussians
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
var_Gaussian <- function(transformation, N, hatd = NULL){
  if(transformation == "cohensd"){
    # sqrt(N)\hat d ~ non-central t with nc = sqrt(N)d, df = N-1
    df = N - 1
    nc = sqrt(N) * hatd

    if( df < 320){
      fac <- sqrt(df / 2) * gamma((df - 1) / 2) / gamma(df / 2)
    }else{
      fac <- 1 / (1 - 3 / (4 * df -1 ))
    }

    v <- ( df * (1 + nc^2) / (df - 2) - nc^2 * fac^2 )
  }
  else if(transformation == "skewness"){
    v <- 6 * (N - 2) / (N + 1) / (N + 3)
  }else if(transformation == "kurtosis"){
    v <- 24 * N * (N - 2) * (N - 3) / (N + 1)^2 / (N + 3) / (N + 5)
  }
  return(v)
}


#' This functions computes the true standard error depending on N for different
#' functionals of Gaussians
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
se_Gaussian <- function(transformation, N, hatd = NULL){
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

#' This functions computes the true standard error depending on N for different
#' functionals of Gaussians
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
gamma1_Gaussian <- function(transformation, N, hatd = NULL){
  if(transformation == "cohensd"){# sqrt(N)\hat d ~ non-central t with nc = sqrt(N)d, df = N-1
    gamma1 <- NaN
  }
  else if(transformation == "skewness"){
    gamma1 <- 0
  }else if(transformation == "kurtosis"){
    gamma1 <- 6 * (N^2 - 5*N + 2) / (N + 7) / (N + 9) *
              sqrt(6 * (N + 3) * (N + 5) / N / (N - 2) / (N - 3))
  }
  return(gamma1)
}

#' This functions computes the true standard error depending on N for different
#' functionals of Gaussians
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
gamma2_Gaussian <- function(transformation, N, hatd = NULL){
  if(transformation == "cohensd"){# sqrt(N)\hat d ~ non-central t with nc = sqrt(N)d, df = N-1
    gamma2 <- NaN
  }
  else if(transformation == "skewness"){
    gamma2 <- 36 * (N - 7) * (N^2 + 2*N - 5) / ((N - 2) * (N + 5) * (N + 7) * (N + 9))
    gamma2 <- 3 * (N^2 +27*N - 70) * (N + 1) * (N + 3) / ((N - 2) * (N + 5) * (N + 7) * (N + 9))
  }else if(transformation == "kurtosis"){
    gamma2 <- NaN
  }
  return(gamma2)
}

#' This functions computes the bias of certain transformations under Gaussianity
#'
#' @param transformation The considered transformation of the moments
#' @param N Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
bias_Gaussian <- function(transformation, N, hatd = NULL){
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
scaleField <- function(Y,
                       x = seq(0, 1, length.out = nrow(Y)),
                       xeval = seq(0, 1, length.out = 2 * nrow(Y)),
                       h = seq(x[2], 0.1, length.out = 20),
                       method = "loclin",
                       kernel = locpol::gaussK, Weights = NULL ){
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
method_gen <- function(name, SCoPEStype, mu1est, N, kN, m0 = NULL, R = NULL){
  arguments = ls()
  q_method <- list(
    name = name,
    SCoPEStype = SCoPEStype,
    mu1Cest    = mu1est,
    kN         = kN
  )
  # q estimation specification
  if(name == "t"){
    q_method$df = N - 1
    q_method$kN = kN
  }else if(name == "mboot"){
    if(!("Boottype" %in% arguments)){
      Boottype = "t"
    }
    if(!("weights" %in% arguments)){
      weights = "rademacher"
    }
    if(!("Mboots" %in% arguments)){
      Mboots = 5e3
    }
    q_method$Boottype   = Boottype
    q_method$weights    = weights
    q_method$Mboots     = Mboots
    q_method$R          = R
    q_method$kN = kN
  }
  if(!is.null(mu1est)){
    if(mu1est == "m0"){
      q_method$m0 = m0
    }
  }

  return(q_method)
}

#' This function generates a list with the required input for different
#' methods to estimate the quantile function q of a SCoPE set.
#'
#' @param name String name of the implemented method to
#' estimate the quantile function
#' @param ... Optional parameters for the quantile estimation method
#' @return q.method
#' @export
q_method_gen <- function(name, ...){
  # Get the input arguments
  arguments = as.list(match.call())

  # Initialize the q_method list, which will be the output
  q.method <- list(
    name = name
  )

  # Fill the q_method list for different quantile function estimation
  # methods
  if(q.method$name %in% c("mboot", "fair.mboot")){
    # Use default values for the multiplier bootstrap, if
    # the values are not provided in the input arguments.
    if(!("Boottype" %in% arguments)){
      Boottype = "t"
    }
    if(!("weights" %in% arguments)){
      weights = "rademacher"
    }
    if(!("Mboots" %in% arguments)){
      Mboots = 5e3
    }

    q.method$Boottype   = Boottype
    q.method$weights    = weights
    q.method$Mboots     = Mboots
    q.method$R          = R

    # Add the fairness parameters if the bootstrap
    # should be fair
    if(q.method$name == "fair.mboot"){
      q.method$fair.intervals = fair.intervals
      if(!("fair.type" %in% arguments)){
        q.method$fair.type = "linear"
      }else{
        q.method$fair.type = fair.type
      }

      if(!("fair.niter" %in% arguments)){
        q.method$fair.niter = 10
      }else{
        q.method$fair.niter = fair.niter
      }
    }
  }else if(q.method$name == "t.iid"){
    q.method$df      = N - 1
  }

  return(q.method)
}

#' This function generates a list with the required input for different
#' methods to estimate the quantile function q of a SCoPE set.
#'
#' @param name String name of the implemented method to
#' estimate the quantile function
#' @param ... Optional parameters for the quantile estimation method
#' @return Scale field
#' @export
Preimage_method_gen <- function(name, ...){
  # Get the input arguments
#  arguments = as.list(match.call())

  # Initialize the Preimage.method list, which will be the output
  Preimage.method <- list(
    name = name
  )

  # Fill the Preimage.method list for different Preimage estimation
  # methods
  if(Preimage.method$name == "thickening"){
    Preimage.method$kN = kN

  }else if(Preimage.method$name == "true"){
    Preimage.method$kN = 0
    Preimage.method$mu = mu
  }else if(Preimage.method$name == "SCB"){
    Preimage.method$kN = 0
  }else if(Preimage.method$name == "Storey.iid"){
    Preimage.method$m0 = m0
  }

  return(Preimage.method)
}


#' @export
plot_col <- function(statistic, C, detect, x = NULL,
                     xlab = '', ylab = '', title = '',
                     mu = NULL){
  if(is.null(x)){
    x = 1:length(statistic)
  }
  y = statistic
  # Get a color vector indicating the "red" (upper excursions) and
  # the blue (lower excursions) set
  colVec <- rep("black", length(y))
  colVec[detect] <- "red"
  plot(x, y, col = colVec,
       pch = 18, xlab = xlab, ylab = ylab, main = title)
  lines(x, C, lty = 2, col = "orchid3")
}
