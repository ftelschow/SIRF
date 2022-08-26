#------------------------------------------------------------------------------#
#                                                                              #
#     Bootstrap functions to obtain the quantiles of the maximum of random     #
#     fields                                                                   #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#
# Contained functions:
#      - Fairthreshold1D
#      - Fair
#
#------------------------------------------------------------------------------#
# Developer notes:
#
#------------------------------------------------------------------------------#

#' Estimates from a sample of random functions (for example a bootstrap sample)
#' a Fair thresholding function q
#'
#' @param sample array of dimension K x N containing N-realizations of
#'  a random field over a 1-dimensional domain.
#' @param x vector of length K of locations at which the sample is observed.
#' @param alpha numeric the targeted upper quantile of the maxiumum of the
#'   absolute value of the random field. Default is 0.95.
#' @param fair a vector partitioning the vector x into regions on which the
#'             on which the rejection should be fair. First element must be
#'             x[1] and last x[length(x)]
#' @return list with elements
#'  \itemize{
#'   \item q Vector containing the fair piecewise linear thresholding function at each x
#'   \item qm Vector containing the offset and the slopes of the fair thresholding function
#'   \item EmpRejections Numeric giving the empirical rejection rate of the fair
#'   thresholding function with respect to the sample.
#' }
#' @export
FairThreshold1D <- function(samples, x, fair, fair.type = "linear", Splus = rep(TRUE, length(x)), Sminus = rep(TRUE, length(x)), alpha = 0.95, diff.fair = NULL, subI = NULL ){
  # Get the length of the different intervals for the fairness
  if(is.null(diff.fair)){
    dfair = diff(fair)
  }else{
    dfair = diff.fair
  }
  # Get the subinterval indices
  if(is.null(subI)){
    subI = list()
    for(k in 2:length(fair)){
      subI[[k-1]] = which( x  >= fair[k-1] & x  <= fair[k] )
    }
  }

  dimS = dim(samples)
  #  if(fair.type = "linear"){
  samples_minus = -samples
  samples_minus[!Sminus,] <- -Inf
  samples_plus = samples
  samples_plus[!Splus,]   <- -Inf
  #}

  # initialize the quantile piecewise linear function
  q0 = -Inf

  # iterate over the different intervals
  for(k in 2:length(fair)){
    # find indices within the k-th interval
    subIk = subI[[ k-1]]

    if(q0 == -Inf){
      # Get the local test statistic
      maxIk = apply(rbind(apply(samples_minus[subIk,], 2, max),
                          apply(samples_plus[subIk,], 2, max)), 2, max)

      if(k == 2){
        # Initialize the slope of the quantile function
        mq = quantile( maxIk, 1 - alpha*dfair[k-1], type = 8 )

        # Initialize the quantile function, which is piecewise linear!
        qx = unlist(rep(mq, length(subIk)))

        # Get back into this loop if type is not linear
        if(fair.type == "constant"){
          q0 = -Inf
        }else{
          q0 = mq
        }

      }else{
        mq <- c( mq, quantile( maxIk, 1 - alpha * dfair[k-1], type = 8 ) )
        qx <- c( qx, rep(mq[k-1], length(subIk)) )

        # Update the quantile function, the slope and the starting point for next iteration
        if(fair.type == "constant"){
          q0 = -Inf
        }else{
          q0 = qx[length(qx)]
        }
      }

    }else{

      maxFun <- matrix(mapply(max, samples_plus[subIk,],
                              samples_minus[subIk,]),
                       length(subIk), dimS[2])

      # define the function, which finds the optimal slope for the correct rejection rate
      solvef <- function(l){
        mean(apply( maxFun - q0 - l *(x[subIk] - x[subIk[1]-1] ) > 0, 2, any)) - alpha * dfair[k-1]
      }
      # optimize the rejection rate function
      if( all(is.infinite(maxFun)) ){
        qk <- list()
        qk$root <- -Inf
      }else{
        qk <- uniroot(solvef, interval = c(-500, 500))
      }
      # Update the quantile function, the slope and the starting point for next iteration
      mq <- c(mq, qk$root)
      qx = c(qx, q0 + qk$root * (x[subIk] - x[subIk[1]-1]))
      q0 = q0 + qk$root * dfair[k-1]
      if(is.nan(q0)){
        q0 = -Inf
      }
    }
  }

  # Interval counter to later fill not important intervals
  if(any(is.infinite(mq))){
    if(fair.type == "linear"){
      qx[is.infinite(qx)] = NA
      if(is.na(qx[1])){
        qx[1] = 0
      }
      Eqx = length(qx)
      if(is.na(qx[Eqx])){
        qx[Eqx] = 0
      }
      qx = na.approx(qx)
    }else{
      qx[is.infinite(qx)] = mean(qx[!is.infinite(qx)])
    }
  }
  # Compute the empirical rejection rate of the quantile function on the given sample
  maxFun <- matrix(mapply(max, samples_plus, samples_minus), dimS[1], dimS[2])
  diffmFunqx = maxFun - qx
  diffmFunqx[which(is.nan(diffmFunqx))] = 0

  EmpRejections = mean(apply(diffmFunqx > 0, 2, any))

  # return the results
  return(list(q = qx, mq = mq, EmpRejections = EmpRejections))
}


#' Estimates from a sample of random functions (for example a bootstrap sample)
#' a Fair thresholding function q
#'
#' @param sample array of dimension K x N containing N-realizations of
#'  a random field over a 1-dimensional domain.
#' @param x vector of length K of locations at which the sample is observed.
#' @param alpha numeric the targeted upper quantile of the maxiumum of the
#'   absolute value of the random field. Default is 0.95.
#' @param fair a vector partitioning the vector x into regions on which the
#'             on which the rejection should be fair. First element must be
#'             x[1] and last x[length(x)]
#' @return list with elements
#'  \itemize{
#'   \item q Vector containing the fair piecewise linear thresholding function at each x
#'   \item qm Vector containing the offset and the slopes of the fair thresholding function
#'   \item EmpRejections Numeric giving the empirical rejection rate of the fair
#'   thresholding function with respect to the sample.
#' }
#' @export
OptimizeFairThreshold1D <- function(samples,
                                    x,
                                    fair,
                                    fair.type = "linear",
                                    Splus = rep(TRUE, length(x)),
                                    Sminus = rep(TRUE, length(x)),
                                    alpha = 0.95,
                                    niter = 10,
                                    subI = NULL,
                                    print.coverage = TRUE ){
  # Get the weighting for the interval, by removing empty intervals
  if(!all(Splus) || !all(Sminus)){
    count = 0
    dd = diff(fair)
    diff.fair = 0*(1:length(dd))
    for(k in 2:length(fair)){
      if(is.null(subI)){
        subIk = which(   fair[k-1] <= x & x  < fair[k] )
        if(k == length(fair)){
          subIk = c(subIk, length(x))
        }
      }else{
        subIk = subI[[k-1]]
      }
      if( any(Splus[subIk]) || any(Sminus[subIk]) ){
        diff.fair[k-1] <- dd[k-1]
        count = count + 1
      }
    }

    diff.fair = diff.fair * (length(diff.fair)) / count
    if(fair.type == "linear"){
      vv = (diff.fair == 0)
      diff.fair[vv] = 0
      if(vv[1] == TRUE){
        diff.fair[vv] = 0
      }
    }
  }else{
    diff.fair = diff(fair)
  }

  # Compute the fair threshold function
  test = FairThreshold1D(samples, x, fair, fair.type = fair.type, Splus = Splus, Sminus = Sminus, alpha, diff.fair = diff.fair, subI = subI )

  # Initialize values for computing the fair threshold function
  count = 0
  breakCond = FALSE
  oldEmp    = test$EmpRejections
  eps = max(alpha * 0.025, 10/dim(samples)[2])
  alpha_new = alpha

  # Loop to remove conservativeness of the fair threshold function
  while(breakCond == FALSE && count < niter){
    diffCoverage = oldEmp - alpha
    if( abs(diffCoverage) > eps ){
      if(count == 0){
        if(diffCoverage < 0){
          a = c(alpha, 3 * alpha)
        }else{
          a = c(0.1 * alpha, alpha)
        }
      }else{
        if(diffCoverage < 0){
          a[1] = alpha_new
          if(a[1] == a[2]){
            a[2] = 2 * a[2]
          }
        }else{
          a[2] = alpha_new
          if(a[1] == a[2]){
            a[1] = 0.1 * a[1]
          }
        }
      }

      # update the global alpha
      alpha_new = mean(a)

      # Get new quantile function
      test   = FairThreshold1D(samples, x, fair, fair.type = fair.type,
                               Splus = Splus, Sminus = Sminus,
                               alpha_new, diff.fair = diff.fair, subI = subI )
      count  = count + 1
      oldEmp = test$EmpRejections
      if(print.coverage){
        print(oldEmp)
      }
    }else{
      breakCond = TRUE
    }
  }

  test$loc_alpha = alpha_new * diff.fair

  # return the results
  return(test)
}
