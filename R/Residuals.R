#------------------------------------------------------------------------------#
#
#     Functions for Residuals
#
#------------------------------------------------------------------------------#
# Contained functions:
#  - DeltaResiduals
#  - get_moments
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline included
#------------------------------------------------------------------------------#

#' Computes functional delta residuals for moment based statistics
#'
#' @param Y array of dimension K_1 x ... x K_D x N containing N-realizations of
#'  a random field over a D-dimensional domain.
#' @param transformation expression of the form f ~ f(mu1, ..., muk ) giving
#' the transformation of the moments. Can be a "linear", "cohensd", "skewness" or
#' "kurtosis". Then standard estimators are used.
#' Current options are "skewness" and "kurtosis".
#' @param moments string of the form c("mu1", ..., "muk") containing the names of
#' the moments contained in the order of transformation input.
#' @return an array representing the delta residuals
#' @export
DeltaMomentResiduals <- function(Y,
                                 transformation,
                                 moments = NULL,
                                 oracle  = NULL ){
  # Get the dimension of the field
  dimY     <- dim(Y)
  N.sample <- dimY[2]
  N.locs   <- dimY[1]

  # Choice of statistic if gradf is a string
  #-----------------------------------------------------------------------------
  lin <- FALSE
  if( is.character(transformation) ){
    if( transformation == "linear" ){
      transformation <- f ~ mu1
      moments <- "mu1"
      lin <- TRUE
    }else if( transformation == "cohensd" ){
      transformation <- f ~ mu1 / sqrt( mu2 - mu1^2 )
      moments <- c("mu1", "mu2")
    } else if( transformation == "skewness" ){# Fisher's g1
      transformation <- f ~ ( mu3 - 3*mu2*mu1 + 2*mu1^3 ) / ( mu2 - mu1^2 )^(3/2)
      moments <- c("mu1", "mu2", "mu3")
    }else if( transformation == "kurtosis" ){# Fisher's g2
      transformation <- f ~ (mu4 - 4 * mu1 * mu3 + 6 * mu1^2 * mu2 - 3 * mu1^4) / (mu2 - mu1^2)^2 - 3
      moments <- c("mu1", "mu2", "mu3", "mu4")
    }else if( transformation == "kurtosis (unbiased)" ){
      transformation <- as.formula(bquote(
                          f ~ (.(N.sample) + 1) * (.(N.sample) - 1) / (.(N.sample) - 2) / (.(N.sample) - 3) *
                            (mu4 - 4 * mu1 * mu3 + 6 * mu1^2 * mu2 - 3 * mu1^4 ) /
                            (mu2 - mu1^2 )^2 -
                            3 * (.(N.sample) - 1)^2 / (.(N.sample) - 2) / (.(N.sample) - 3)
                        ))
      moments <- c("mu1", "mu2", "mu3", "mu4")
    }else if( transformation == "skewness (normality)"){
      #---- Skewness transformation
      m21     <- var_Gaussian("skewness", N.sample)
      gamma21 <- gamma2_Gaussian("skewness", N.sample)

      W2    <- sqrt(2 * gamma21 + 4) - 1
      delta1 <- 1 / sqrt(log(sqrt(W2))) / sqrt(N.sample)
      alpha <- sqrt(2 / (W2 - 1))

      transformation <- as.formula(bquote(f ~ .(delta1) * log( ( mu3 - 3*mu2*mu1 + 2*mu1^3 ) /
                                                      ( mu2 - mu1^2 )^(3/2) /
                                                      sqrt(.(m21)) / .(alpha) +
                                                      sqrt(  ( mu3 - 3*mu2*mu1 + 2*mu1^3 )^2 /
                                                               ( mu2 - mu1^2 )^(3) /
                                                               .(m21) / .(alpha^2) + 1) )))

      moments <- c("mu1", "mu2", "mu3")
    }else if(transformation == "kurtosis (normality)"){
      # kurtosis transformation
      m12 <- mean_Gaussian("kurtosis", N.sample)
      m22 <- var_Gaussian("kurtosis", N.sample)
      gamma12 <- gamma1_Gaussian("kurtosis", N.sample)

      A <- 6 + 8 / gamma12 * (2 / gamma12 + sqrt(1 + 4 / gamma12^2))
      delta2 <- sqrt(9 * A / 2 / N.sample)

      transformation <- as.formula(bquote(f ~  .(delta2) * ((1 - 2 / 9 / .(A)) - ( (1 - 2 / .(A)) / (1 + ((mu4 - 4 * mu1 * mu3 + 6 * mu1^2 * mu2 - 3 * mu1^4) / (mu2 - mu1^2)^2 - .(m12))  / sqrt(.(m22)) * sqrt(2 / (.(A) - 4))) )^(1/3) )))

      moments = c("mu1", "mu2", "mu3", "mu4")
    }else if(transformation == "Ksquare"){
      #---- Skewness transformation
      m21     <- var_Gaussian("skewness", N.sample)
      gamma21 <- gamma2_Gaussian("skewness", N.sample)
      W2    <- sqrt(2 * gamma21 + 4) - 1
      delta1 <- 1 / sqrt(log(sqrt(W2))) / sqrt(N.sample)
      alpha <- sqrt(2 / (W2 - 1))

      # kurtosis transformation
      m12 <- mean_Gaussian("kurtosis", N.sample)
      m22 <- var_Gaussian("kurtosis", N.sample)
      gamma12 <- gamma1_Gaussian("kurtosis", N.sample)
      A <- 6 + 8 / gamma12 * (2 / gamma12 + sqrt(1 + 4 / gamma12^2))
      delta2 <- sqrt(9 * A / 2 / N.sample)

      transformation <-  as.formula(bquote(f ~  .(delta2^2) * ((1 - 2 / 9 / .(A)) - ( (1 - 2 / .(A)) /
                                                                             (1 + ((mu4 - 4 * mu1 * mu3 + 6 * mu1^2 * mu2 - 3 * mu1^4) /
                                                                                     (mu2 - mu1^2)^2 - .(m12))  / sqrt(.(m22)) *
                                                                                sqrt(2 / (.(A) - 4))) )^(1/3) )^2 +
                                  .(delta1^2) * log( ( mu3 - 3*mu2*mu1 + 2*mu1^3 ) /
                                                       ( mu2 - mu1^2 )^(3/2) /
                                                       sqrt(.(m21)) / .(alpha) +
                                                       sqrt(  ( mu3 - 3*mu2*mu1 + 2*mu1^3 )^2 /
                                                                ( mu2 - mu1^2 )^(3) /
                                                                .(m21) / .(alpha^2) + 1) )^2 ))
      moments = c("mu1", "mu2", "mu3", "mu4")
    }else if(transformation == "Ksquare (normality)"){
      #---- Skewness transformation
      m21     <- var_Gaussian("skewness", N.sample)
      gamma21 <- gamma2_Gaussian("skewness", N.sample)
      W2      <- sqrt(2 * gamma21 + 4) - 1
      delta1  <- 1 / sqrt(log(sqrt(W2)))
      alpha   <- sqrt(2 / (W2 - 1))

      # kurtosis transformation
      m12 <- mean_Gaussian("kurtosis", N.sample)
      m22 <- var_Gaussian("kurtosis", N.sample)
      gamma12 <- gamma1_Gaussian("kurtosis", N.sample)
      A <- 6 + 8 / gamma12 * (2 / gamma12 + sqrt(1 + 4 / gamma12^2))
      delta2 <- sqrt(9 * A / 2)
      p = 3
      transformation <-  as.formula(bquote(f ~  (((.(delta2^2) * ((1 - 2 / 9 / .(A)) - ( (1 - 2 / .(A)) /
                                                                                        (1 + ((mu4 - 4 * mu1 * mu3 + 6 * mu1^2 * mu2 - 3 * mu1^4) /
                                                                                                (mu2 - mu1^2)^2 - .(m12))  / sqrt(.(m22)) *
                                                                                           sqrt(2 / (.(A) - 4))) )^(1/3) )^2 / 2 +
                                             .(delta1^2) * log( ( mu3 - 3*mu2*mu1 + 2*mu1^3 ) /
                                              ( mu2 - mu1^2 )^(3/2) /
                                              sqrt(.(m21)) / .(alpha) +
                                              sqrt(  ( mu3 - 3*mu2*mu1 + 2*mu1^3 )^2 /
                                              ( mu2 - mu1^2 )^(3) /
                                    .(m21) / .(alpha^2) + 1) )^2 / 2)^(1/.(p)) - (1 - 1 / .(p^2) )) / sqrt(.(N.sample) / .(p^2)) - 0) /  1 ))#0.01233952) /  0.9736893 ))
      moments = c("mu1", "mu2", "mu3", "mu4")
    }
  }

  # Get the number of moments, i.e., the dimension of the domain of the
  # transformation
  N.moments <- length(moments)

  #----- Get the untransformed residuals
  #-----------------------------------------------------------------------------
  # Get the moment residuals and all its derivates
  tmp = get_moments(Y, moments)
  # Estimates of the moments
  hatmu  = tmp$moments

  if( is.matrix(oracle) ){
    mu = oracle
  }else{
    mu = hatmu
  }

  # Moment residuals
  res = tmp$moments.residuals

  # Get the estimates of covariance matrix of moments to correct for bias
  hatcov.moment  <- tmp$hatcov.moment
  hess.weights   <- tmp$hess.weights
  N.hess.weights <- length(hess.weights)

  # Get the mixed residuals for second order approximations
  mixed.moment.residuals <- tmp$mixed.moment.residuals
  rm( tmp )

  #----- compute the transformation dependend parts of the Taylor expansion
  #-----------------------------------------------------------------------------
  # Get the symbolic derivates up to order 2
  gradf <- deriv( transformation, moments, func = TRUE, hessian = TRUE )

  # value of the statistic
  stat <- apply( mu, 1, function(x) do.call(gradf, as.list(x)) )

  # gradient of the transformation evaluated at the estimated moments
  gradf.mu <- t(apply( mu,
                        1,
                        function(x) attr(do.call(gradf, as.list(x)), "gradient") ))

  if(lin){
    gradf.mu <- t(gradf.mu)
  }

  # hessian of the transformation evaluated at the estimated moments
  hessf.mu <- matrix( NaN, N.locs, N.hess.weights )
  for( n.l in 1:N.locs ){
    HessM <- matrix( attr(do.call(gradf, as.list(mu[n.l, ])), "hessian"),
                     N.moments,
                     N.moments )
    hessf.mu[n.l, ] <- HessM[ lower.tri( HessM, diag = T ) ] * hess.weights
  }
  rm(n.l, HessM)

  #----- compute the delta residuals
  # First order delta residuals
  delta.res <- t(vapply( 1:N.locs,
                         FUN = function(l) gradf.mu[l, ] %*% res[l, , ],
                         FUN.VALUE = rep(NaN, N.sample) ))

  # Add the second order term if necessary
  delta.res.order2 <- 1 / (2*sqrt(N.sample)) *
                        vapply( 1:N.sample,
                                  FUN = function(n) rowSums( hessf.mu * mixed.moment.residuals[, , n] ),
                                  FUN.VALUE = rep(NaN, N.locs) )

  #----- simple estimate the bias
  hatbias = rowSums( hessf.mu * hatcov.moment ) / 2

  #----- Get derived quantities
  # get the standard deviation of the residuals
  sd.delta.res <- apply( delta.res, 1, sd)
  # get the standard deviation of the residuals
  delta.res.order2 <- delta.res - delta.res.order2
  sd.delta.res.order2 <- apply( delta.res.order2, 1, sd)

  #----- output the delta results and the standard deviation
  return( list( delta.res        = delta.res,
                delta.res.order2 = delta.res.order2,
                delta.sd         = sd.delta.res,
                delta.sd.order2  = sd.delta.res.order2,
                hatmoments       = hatmu,
                moments.res      = res,
                statistic         = stat,
                hatbias           = hatbias,
                gradf.hatmoments  = gradf.mu,
                hessf.hatmoments  = hessf.mu ) )
}


#' Computes the non-central moments of a functional data sets
#'
#' @param Y matrix containing the data. Last dimension enumerates the samples.
#' @param moments string of the form c("mu1", ..., "muk") containing the names.
#' @param cmoments boolean indicating whether an estimate of the cvariance matrix
#' of the moment vector is included in the output
#' @param squeeze boolean
#' @return an list containing the moments and the corresponding moment residuals.
#' #' @export
get_moments <- function( Y, moments ){
  # Get the dimension of the field
  dimY     <- dim(Y)
  N.sample <- dimY[2]
  N.locs   <- dimY[1]

  # Make moments input vector into numeric vector
  moments <- vapply( moments,
                     function(x) as.numeric(substr(x, 3, 3)),
                     FUN.VALUE = 1 )
  N.moments <- length(moments)

  #----- Compute moments and moment residuals
  mu  <- matrix(NaN, N.locs, N.moments)
  res <- array(NaN, dim = c(N.locs, N.moments, N.sample))
  for( k in 1:N.moments ){
    Yk       <- Y^moments[k]
    mu[, k]  <- rowMeans(Yk)
    res[, k,] <- Yk - mu[,k]
  }
  rm(k, Yk)

  # Get the needed cross product terms
  cvec = NULL
  for( k in 1:N.moments ){
    cvec <- rbind( cvec, cbind( k, k:N.moments ) )
  }
  rm(k)
  N.cvec <- dim(cvec)[1]

  #----- Compute estimate of covariance of moments and the mixed moment residuals
  hatcov.moment   <- array(NaN, dim = c(N.locs, N.cvec) )
  mixed.residuals <- array(NaN, dim = c(N.locs, N.cvec, N.sample) )
  hess.weights    <- rep(2, N.cvec)
  for( k in 1:N.cvec ){
    hatcov.moment[,k]  <- rowMeans( Y^(moments[cvec[k, 1]] + moments[cvec[k, 2]]) ) -
                                mu[, cvec[k, 1]] * mu[, cvec[k, 2]]
    # Change weights to one, if this element is unique in the covariance matrix
    if( cvec[ k, 1 ] == cvec[ k, 2 ] ){
      hess.weights[k] = 1
    }

    mixed.residuals[,k,] <- res[, cvec[k, 1], ] * res[, cvec[k, 2], ]
  }

  #----- return list containing the moments and the moment residuals
    return( list(moments           = mu,
                 hatcov.moment     = hatcov.moment / N.sample,
                 hess.weights      = hess.weights,
                 hess.entries      = cvec,
                 moments.residuals = res,
                 mixed.moment.residuals = mixed.residuals) )
}
