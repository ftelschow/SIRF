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
FairThreshold1D <- function(samples, x, fair.intervals,
                            fair.type = "linear",
                            crit.set,
                            alpha = 0.05,
                            diff.fair = NULL,
                            subI = NULL, inter = NULL ){
  #
  if(is.null(subI) || is.null(inter)){
    s = sub.intervals(x, fair.intervals, crit.set)
    subI  = s$subI
    inter = s$inter
  }

  #
  Sminus = crit.set$minus
  Splus  = crit.set$plus

  # Get the length of the different intervals for the fairness
  if(is.null(diff.fair)){
    dfair = diff(fair.intervals)
  }else{
    dfair = diff.fair
  }

  # Spread the parts where no probability is expected fairly
  # to the other partitions
  if(!all(inter)){
    dfair = dfair + sum(dfair[!inter]) / sum(inter)
    dfair[!inter] = 0
  }

  dimS = dim(samples)
  #  if(fair.type = "linear"){
  samples_minus = -samples
  samples_minus[!Sminus,] <- -Inf
  samples_plus = samples
  samples_plus[!Splus,]   <- -Inf
  #}

  subI_prob <- function(k, q){
    if(is.null(crit.set$minus)){
      low.excursion = rep(FALSE, dim(samples)[2])
    }else{
      ind_minus     = intersect(which(crit.set$minus), subI[[k]])
      if(length(ind_minus) == 0){
        low.excursion = rep(FALSE, dim(samples)[2])
      }else if(length(ind_minus) == 1){
        low.excursion = -samples[ind_minus,] - q[ind_minus] >= 0
      }else{
        low.excursion = apply(-samples[ind_minus,] - q[ind_minus], 2, max) >= 0
      }
    }

    if(is.null(crit.set$plus)){
      up.excursion = rep(FALSE, dim(samples)[2])
    }else{
      ind_plus     = intersect(which(crit.set$plus), subI[[k]])
      if(length(ind_plus) == 0){
        up.excursion = rep(FALSE, dim(samples)[2])
      }else if(length(ind_plus) == 1){
        up.excursion = samples[ind_plus,] - q[ind_plus] >= 0
      }else{
        up.excursion = apply(samples[ind_plus,] - q[ind_plus], 2, max) >= 0      }
    }

    mean(apply(cbind(low.excursion, up.excursion), 1, any))
  }

  # initialize the quantile piecewise linear function
  q0 = -Inf
  q  = rep(-Inf, length(x))

  # Get the indices of the partitions which are used
  # to control the limit distribution
  Ik = c(which(inter), -Inf)

  # iterate over the different intervals
  for(k in 1:length(Ik[-length(Ik)])){
    # find indices within the k-th interval
    subIk = subI[[Ik[k]]]

    if(q0 == -Inf){
      # Get the local test statistic
      maxIk = apply(rbind(apply(samples_minus[subIk,], 2, max),
                          apply(samples_plus[subIk,], 2, max)), 2, max)

      # Initialize the piecewise linear function
      mq = quantile( maxIk, 1-alpha*dfair[Ik[k]], type = 8 )

      q[subIk] = mq

      # Get back into this loop if type is not linear
      if(fair.type == "constant" || Ik[k]+1 != Ik[k+1]){
        q0 = -Inf
      }else{
        q0 = mq
      }
    }else{
      # define the function, which finds the optimal slope for the correct rejection rate
      solvef <- function(l){
        qq = q
        qq[subIk] = q0 - l * (x[subIk] - x[subIk[1]-1])

        return(subI_prob(Ik[k], qq) - alpha * dfair[Ik[k]])
      }
      qk <- uniroot(solvef, interval = c(-500, 500))

      q[subIk] = q0 - qk$root * (x[subIk] - x[subIk[1]-1])

      # Get back into this loop if type is not linear
      if(fair.type == "constant" || Ik[k]+1 != Ik[k+1]){
        q0 = -Inf
      }else{
        q0 = q[subIk[length(subIk)]]
      }
    }
  }

  # Interval counter to later fill not important intervals
  if(any(is.infinite(q))){
    if(fair.type == "linear"){
      q[is.infinite(q)] = NA
      if(is.na(q[1])){
        q[1] = mean(q, na.rm = TRUE)
      }
      Eqx = length(q)
      if(is.na(q[Eqx])){
        q[Eqx] = mean(q, na.rm = TRUE)
      }
      q = na.approx(q)
    }else{
      q[is.infinite(q)] = mean(q[!is.infinite(q)])
    }
  }

  EmpRejections = IntervalProb(q, crit.set, samples, x, fair.intervals, subI = subI)

  # return the results
  return(list(q = q, mq = mq, EmpRejections = EmpRejections))
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
                                    fair.intervals,
                                    fair.type = "linear",
                                    crit.set,
                                    alpha  = 0.05,
                                    niter  = 10,
                                    diff.fair = NULL,
                                    subI   = NULL,
                                    inter = NULL,
                                    print.coverage = TRUE ){
  # Fill subI and inter, if not provided
  if(is.null(subI) || is.null(inter)){
    s = sub.intervals(x, fair.intervals, crit.set)
    subI  = s$subI
    inter = s$inter
  }

  # Get the stepsize for the interval iteration
  if( sum(crit.set$minus)/2 + sum(crit.set$minus)/2 > 0.75*length(x)  ){
    stepsize = 3
  }else{
    stepsize = 1.5
  }

  # Compute the fair threshold function
  test = FairThreshold1D(samples = samples,
                         x = x,
                         fair.intervals = fair.intervals,
                         fair.type = fair.type,
                         crit.set = crit.set,
                         alpha = alpha,
                         diff.fair = diff.fair,
                         subI = subI, inter = inter )

  # Initialize values for computing the fair threshold function
  count = 0
  breakCond = FALSE
  oldEmp    = test$EmpRejections$global
  eps = max(alpha * 0.05, 10/dim(samples)[2])
  alpha_new = alpha

  # Loop to remove conservativeness of the fair threshold function
  while(breakCond == FALSE && count < niter){
    diffCoverage = oldEmp - alpha
    if( abs(diffCoverage) > eps ){
      if(count == 0){
        if(diffCoverage < 0){
          a = c(alpha, stepsize * alpha)
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
      test   = FairThreshold1D(samples = samples,
                               x = x,
                               fair.intervals = fair.intervals,
                               fair.type = fair.type,
                               crit.set = crit.set,
                               alpha = alpha_new,
                               diff.fair = diff.fair,
                               subI = subI, inter = inter )
      count  = count + 1
      oldEmp = test$EmpRejections$global
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


#' Estimates using a sample of random functions (for example a bootstrap sample)
#' the FWER on different intervals.
#'
#' @param sample array of dimension K x N containing N-realizations of
#'  a random field over a 1-dimensional domain.
#' @param x vector of length K of locations at which the sample is observed.
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
IntervalProb <- function(q, crit.set, samples, x, fair.intervals, subI = NULL){
  # Get the subinterval indices
  if(is.null(subI)){
    subI = list()
    for(k in 2:length(fair.intervals)){
      subI[[k-1]] = which(x  >= fair.intervals[k-1] & x  <= fair.intervals[k])
    }
  }

  subI_prob <- function(k){
    if(is.null(crit.set$minus)){
      low.excursion = rep(FALSE, dim(samples)[2])
    }else{
      ind_minus     = intersect(which(crit.set$minus), subI[[k]])
      if(length(ind_minus) == 0){
        low.excursion = rep(FALSE, dim(samples)[2])
      }else if(length(ind_minus) == 1){
        low.excursion = -samples[ind_minus,] - q[ind_minus] >= 0
      }else{
        low.excursion = apply(-samples[ind_minus,] - q[ind_minus], 2, max) >= 0
      }
    }

    if(is.null(crit.set$plus)){
      up.excursion = rep(FALSE, dim(samples)[2])
    }else{
      ind_plus     = intersect(which(crit.set$plus), subI[[k]])
      if(length(ind_plus) == 0){
        up.excursion = rep(FALSE, dim(samples)[2])
      }else if(length(ind_plus) == 1){
        up.excursion = samples[ind_plus,] - q[ind_plus] >= 0
      }else{
        up.excursion = apply(samples[ind_plus,] - q[ind_plus], 2, max) >= 0
      }
    }

    mean(apply(cbind(low.excursion, up.excursion), 1, any))
  }

  # Intervalwise rejections
  subI.probs = vapply(1:length(subI), subI_prob, FUN.VAL = 0.1)

  # Global rejectionrate
  if(is.null(crit.set$minus)){
    low.excursion = rep(FALSE, dim(samples)[2])
  }else{
    low.excursion = apply(-samples[crit.set$minus,] - q, 2, max) >= 0
  }
  if(is.null(crit.set$plus)){
    up.excursion = rep(FALSE, dim(samples)[2])
  }else{
    up.excursion = apply(samples[crit.set$plus,] - q, 2, max) >= 0
  }
  global.prob = mean(apply(cbind(low.excursion, up.excursion), 1, any))

  return(list(global = global.prob, local = subI.probs))
}

#' Estimates using a sample of random functions (for example a bootstrap sample)
#' the FWER on different intervals.
#'
#' @param sample array of dimension K x N containing N-realizations of
#'  a random field over a 1-dimensional domain.
#' @param x vector of length K of locations at which the sample is observed.
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
sub.intervals <- function(x, fair.intervals, crit.set){
  # Get the subintervals
  subI = list()
  for(k in 2:length(fair.intervals)){
    subI[[k-1]] = which(x >= fair.intervals[k-1] & x <= fair.intervals[k])
  }

  # Get the intervals which have an intersection with the critical sets
  subI.cit <- unlist(lapply(subI, function(l){
    ifelse(length(intersect(l, which(crit.set$minus))) != 0 |
             length(intersect(l, which(crit.set$plus))) != 0, TRUE, FALSE )

  } ))

  list(subI = subI, inter = subI.cit)
}
