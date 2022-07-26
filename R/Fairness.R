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
FairThreshold1D <- function(samples, x, fair, fair.type = "linear", Splus = rep(TRUE, length(x)), Sminus = rep(TRUE, length(x)), alpha = 0.95 ){
  # Get the length of the different intervals for the fairness
  dfair = diff(fair)
  dimS = dim(samples)

    # initialize the quantile piecewise linear function
    q0 = -Inf
    # iterate over the different intervals
    for(k in 2:length(fair)){
      if(q0 == -Inf){
        # Get the indices of x where the first interval is contained
        subIk  = which( x  >= fair[k-1] & x  <= fair[k] )
        xx = x
        xx[!Splus] = Inf
        subIk_Splus  = which( xx  >= fair[k-1] & xx  <= fair[k] )
        xx = x
        xx[!Sminus] = Inf
        subIk_Sminus = which( xx >= fair[k-1] & xx <= fair[k] )

        # Get the maximum over the first interval
        if(length(subIk_Splus) == 1){
          maxIk_Splus <- samples[subIk_Splus,]
        }else if(length(subIk_Splus) > 1){
          maxIk_Splus <- apply( samples[subIk_Splus,], 2, max )
        }else{
          maxIk_Splus <- rep(-Inf, Mboots)
        }

        if(length(subIk_Sminus) == 1){
          maxIk_Sminus <- -samples[subIk_Sminus,]
        }else if(length(subIk_Sminus) > 1){
          maxIk_Sminus <- apply( -samples[subIk_Sminus,], 2, max )
        }else{
          maxIk_Sminus <- rep(-Inf, Mboots)
        }

        maxIk  <- apply( rbind(maxIk_Splus, maxIk_Sminus), 2, max )

        if(k == 2){
          # Initialize the slope of the quantile function
          mq = quantile( maxIk, 1 - alpha*dfair[1], type = 8 )

          # Initialize the quantile function, which is piecewise linear!
          qx = unlist(rep( mq, length(subIk)))

          # Get back into this loop if type is not linear
          if(fair.type == "constant"){
            q0 = -Inf
          }else{
            q0 = mq
          }

        }else{
          mq <- c( mq, quantile( maxIk, 1 - alpha*dfair[k-1], type = 8 ) )
          qx <- c( qx, rep(mq[k-1], length(subIk)) )

          # Update the quantile function, the slope and the starting point for next iteration
          if(fair.type == "constant"){
            q0 = -Inf
          }
        }

      }else{
        # find indices within the k-th interval
        subIk  = which( x  >= fair[k-1] & x  <= fair[k] )
        #samples_Ik = samples[subIk,]

        xx = x
        xx[!Splus] = Inf
        PlusSet <- ( xx  >= fair[k-1] & xx  <= fair[k] )
        subIk_Splus  = which( PlusSet )
        xx = x
        xx[!Sminus] = Inf
        MinusSet <- ( xx >= fair[k-1] & xx <= fair[k] )
        subIk_Sminus = which( MinusSet )

          sample_minus = -samples
          sample_minus[!MinusSet,] <- -Inf
          sample_plus = samples
          sample_plus[!PlusSet,]   <- -Inf
          maxFun <- matrix(mapply(max, sample_plus, sample_minus), dimS[1], dimS[2])
          maxFun = maxFun[subIk,]

          # define the function, which finds the optimal slope for the correct rejection rate
          solvef <- function(l){
            mean(apply( maxFun - q0 - l *(x[subIk] - x[subIk[1]-1] ) > 0, 2, any)) - alpha * dfair[k-1]
          }
          # optimize the rejection rate function
          if( all(!PlusSet) && all(!MinusSet) ){
            qk <- list()
            qk$root <- -Inf
          }else{
            qk <- uniroot(solvef,interval = c(-100, 100))
          }
          # Update the quantile function, the slope and the starting point for next iteration
          mq <- c(mq, qk$root)
          qx = c(qx, q0 + qk$root * (x[subIk] - x[subIk[1]-1]))
          q0 = q0 + qk$root * dfair[k-1]

      }
    }

    # Compute the empirical rejection rate of the quantile function on the given sample
    sample_minus = -samples
    sample_minus[!Sminus,] <- -Inf
    sample_plus = samples
    sample_plus[!Splus,]   <- -Inf
    maxFun <- matrix(mapply(max, sample_plus, sample_minus), dimS[1], dimS[2])
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
OptimizeFairThreshold1D <- function(samples, x, fair, fair.type = "linear", Splus = rep(TRUE, length(x)), Sminus = rep(TRUE, length(x)), alpha = 0.95, niter = 10, print.coverage = TRUE ){
  # Compute the fair threshold function
  test = FairThreshold1D(samples, x, fair, fair.type = fair.type, Splus = Splus, Sminus = Sminus, alpha )

  # Initialize values for computing the fair threshold function
  count = 0
  breakCond = FALSE
  oldEmp    = test$EmpRejections

  # Loop to remove conservativeness of the fair threshold function
  while(breakCond == FALSE && count <= niter){
    diffCoverage = oldEmp - alpha
    if(  abs(diffCoverage) > alpha * 0.025 ){
      if(count == 0){
        if(diffCoverage < 0){
          a = c(alpha, 3 * alpha)
        }else{
          a = c(0.5 * alpha, alpha)
        }
      }else{
        if(diffCoverage < 0){
          a[1] = alpha_new
        }else{
          a[2] = alpha_new
        }
      }

      # update the global alpha
      alpha_new = mean(a)

      # Get new quantile function
      test   = FairThreshold1D(samples, x, fair, fair.type = fair.type, Splus = Splus, Sminus = Sminus, alpha_new )
      count  = count + 1
      oldEmp = test$EmpRejections
      if(print.coverage){
        print(oldEmp)
      }
    }else{
      breakCond = TRUE
    }
  }
  # return the results
  return(test)
}
