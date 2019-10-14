################################################################################################
####                                                                                        ####
####             Functions to simulate the covering rates of different SCBs                 ####
####                                                                                        ####
################################################################################################
#
# Included functions:
#   - covRate_1sample (tested)
#   - covRate_2sample (tested)
#
################################################################################################
#' Simulates the covering rate of different methods for constructiong of simultaneous confidence bands (SCBs) for the mean of functional signal plus noise model (SpN) with different error processes.
#'
#' @param trials Integer number of trials to estimate the covering rate.
#' @param scenario String specifying the scenario. Options are:
#' \itemize{
#'   \item SpN:   Signal plus Noise Model with one or two dimensional domain. Make sure that the functions 'mu', 'sigma' and 'noise' produce a correct and consistent output for the dimension you want to simulate.
#'   \item Scale1d: Scale field of Signal plus Noise Model with one dimensional domain
#' }
#' @param param_scenario list containing parameters for the scenario.
#' \itemize{
#'   \item SpN:     'SmoothWeights': Either NULL or array of size Keval x length(x).
#'   \item Scale1d: Scale field of Signal plus Noise Model with one dimensional domain
#' }
#' @param method String the method to compute the SCBs. Options are "tGKF", "GKF", "Bootstrap", "Bootstrapt", "MultiplierBootstrap", "MultiplierBootstrapt". Default value is "tGKF".
#' @param level Numeric the targeted covering probability. Must be strictly between 0 and 1. Default is 0.95.
#' @param N Integer number of samples used to compute the SCBs.
#' @param mu Function computing the population mean of the SpN.
#' @param noise Function generating random N realizations of a noise process. The output must be an array of dimension K1 x ... x Kd x N.
#' @param sigma Function defining the pointwise variance in the SpN.
#' @param sd_ObsNoise Numeric standard deviation of zero mean Gaussian observation noise, which will be added to the SpNM. Default value is 0.
#' @param x Vector of length K containing the grid points, on which the simulated processes are evaluated.
#' @param ... Additional parameter for 'FunctionalDataSample()', i.e. parameter for the noise process.
#' @return list with elements
#'  \itemize{
#'   \item cov.rate Number (simulated covering rate)
#'   \item trials
#'   \item BootsReplicats
#'   \item nSamp
#' }
#' @export
covRate_simulation <- function(
  trials, scenario, param_scenario=NULL, method, param_method=NULL, level=0.95, N, mu=function(x){ rep(0, length(x)) }, noise=GaussDensitySumNoise, sigma=function(x){rep(1, length(x))}, sd_ObsNoise=0, x=seq(0,1,length.out=100), ...
){
  ######## Check input
  ## trials
  if( is.numeric(trials) ){
    if( trials%%1 !=0 & trials<=0 ){
      stop("The input 'trials' needs to be a positiv natural number.")
    }
  }else{ stop("The input 'trials' needs to be a positiv natural number.") }

  ## Check input scenario
  if( is.character(scenario) ){
    if( !(scenario%in%c("SpN", "Scale1d", "SNR")) ){
      stop("Choose a valid 'scenario' for the simulation. Options are 'SpN', 'Scale1d', 'SNR'.")
    }
  }else{
    stop("The input 'scenario' must be a string. Please, check the help page for available options.")
  }

  ##### Check or initialize the param_method list depending whether it s provided or not.
  ## Signal plus noise model case
  if( scenario=="SpN" ){
    if( is.null(param_scenario) ){
      param_scenario      = list()
      # Put the standard value for the bootstrap replicates
      param_scenario$SmoothWeights = NULL
    }else{
      if( !is.matrix(param_scenario$SmoothWeights) ){
        stop("The input 'SmoothWeights' in the 'param_scenario' list needs to be a matrix.")
      }
    }
  }

  ## Scale space case
  if( scenario=="Scale1d" ){
    if( is.null(param_scenario) ){
      param_scenario        = list()
      # Put the standard value for the LKC estimator
      param_scenario$bw     = LKC_estim_direct
      param_scenario$evalN  = 200
    }else{
      if( is.function(param_scenario$LKC_estim) ){
        stop("The input 'LKC_estim' must be a function computing the LKCs from a sample of random fields.")
      }
    }
  }

  # Check input method
  # if( !is.character(method) ){
  #   stop("The input 'method' must be a string. Please, check the help page for available options.")
  # }else{if( !(method%in%c("tGKF", "GKF", "NonParametricBootstrap", "Bootstrapt", "MultiplierBootstrap", "MultiplierBootstrapt")) ){
  #   stop("Choose a valid option from the available quantile approximations. Please, check the help page.")
  # }}
  # Check input level
  if( !is.numeric(level) ){
    stop("The input 'level' needs to be strictly between 0 and 1.")
  }else{if( level>1 | level<0 ){
    stop("The input 'level' needs to be strictly between 0 and 1.")
  }}
  # Check input N
  # if( !is.numeric(N) ){
  #   stop("The input 'N' needs to be a positiv natural number.")
  # }else{if( N%%1 !=0 & N<=0 ){
  #   stop("The input 'N' needs to be a positiv natural number.")
  # }}
  # Check input x
  if( !is.vector(x) ){
    stop("The input 'x' needs to be a vector.")
  }
  # Check input mu
  if( !is.function(mu) ){
    stop("The input 'mu' needs to be a function.")
  }
  # Check input noise
  if( !is.function(noise) ){
    stop("The input 'noise' needs to be a function.")
  }
  # Check input sigma
  if( !is.function(sigma) ){
    stop("The input 'sigma' needs to be a function.")
  }
  # Check input sd_obsNoise
  if( !is.numeric(sd_ObsNoise) ){
    stop("The input 'sd_ObsNoise' needs to be a positiv real number or zero.")
  }else{if( !is.numeric(sd_ObsNoise) | sd_ObsNoise<0 ){
    stop("The input 'sd_ObsNoise' needs to be a positiv real number or zero.")
  }}

  ###################### Simulation of the covering rate
  # Check whether smoothing is to be applied
  smoothTrue = !is.null(param_scenario$SmoothWeights)

  ##### SpN scenario
  if( scenario %in%  c("SpN", "SNR")  ){
    ##### Compute the mean on the sampling grid for checking whether the SCBs cover it, and initialize the covering rates
    rate           = matrix(rep(0, length(method)*length(N) ),  length(method), length(N))
    colnames(rate) = N
    rownames(rate) = method
    K              = length(x)
    if( scenario == "SpN"){
      mu.check = mu(x)
    }else{
      mu.check = mu(x) / sigma(x)
    }

    if(is.vector(mu.check)){
      D = 1
    }else{
      D = length(dim(mu.check))
    }
    ##### True model parameters for smoothed model
    if(smoothTrue){
      K                = dim(param_scenario$SmoothWeights)[1]
      mu.check         = mu(x)
      musmoothed.check = as.vector(param_scenario$SmoothWeights%*%mu.check)
      mu.check         = mu(param_scenario$xeval)
      rate_smoothed    = matrix(rep(0, length(method)*length(N) ),  length(method), length(N))
      colnames(rate_smoothed) = N
      rownames(rate_smoothed) = method
    }

    ### loop over simulations
    for( i in 1:trials ){
      ### Generate random samples
      Y   = FunctionalDataSample( N=max(N), x=x, mu=mu, noise=noise, sigma=sigma, sd_ObsNoise=sd_ObsNoise)#,... )

      if( smoothTrue ){
        Y = param_scenario$SmoothWeights%*%Y
      }

      ###### Loop over the different methods
      for(n in 1:length(N)){
        if( D==1 ){
          Ytmp = Y[,1:N[n]]
        }else{
          Ytmp = Y[,,1:N[n]]
        }

        for( count_method in 1:length(method) ){
          ### Compute the SCBs
          if( scenario == "SpN"){
            scb = scb_mean( Ytmp, level=level, method=method[count_method], param_method=param_method[[count_method]] )$scb
          }else{
            scb = scb_SNR( Ytmp, level=level, method=method[count_method], param_method=param_method[[count_method]] )$scb
          }

          ### Check whether the SCB covers the true mean or respective the smoothed mean everywhere
          if( all(scb$lo<=mu.check) & all(scb$up>=mu.check) ){
            rate[count_method,n] <- rate[count_method,n] + 1
          }
          if(!is.null(param_scenario$SmoothWeights)){
            if( all(scb$lo<=musmoothed.check) & all(scb$up>=musmoothed.check) ){
              rate_smoothed[count_method,n] <- rate_smoothed[count_method,n] + 1
            }
          }

       }
      }
    }
  }

  ##### Scale1d scenario
  if( scenario == "Scale1d"){
    ##### Get or initialize values for simulation
    K       = length(param_scenario$xeval)
    mu.sim  = mu(x)
    H       = length(param_scenario$bw)

    ##### Initialize the smoothing weights, the smoothed mean and the covering rate counter
    Weights  = array( NA, c(K, length(x), H) )
    mu.check = array( NA, c(K, H) )
    rate     = 0

    ##### Compute the smoothing array for creating the scale field
    for(j in 1:H){
      Weights[,,j]  <- locpol::locLinWeightsC( x=x, xeval=param_scenario$xeval, bw=param_scenario$bw[j], kernel=param_scenario$kernel )$locWeig
    }
    ##### Compute the scale field version of the population mean
    mu.check = scaleField( mu.sim, Weights=Weights )[,,1]

    ##### loop over simulations
    for(i in 1:trials){
      ### Generate random samples and transform them into a scale field
      Y   = FunctionalDataSample( N=N, x=x, mu=mu, noise=noise, sigma=sigma, sd_ObsNoise=sd_ObsNoise,... )
      Y   = scaleField( Y, Weights=Weights )
      ### Compute the SCBs
      scb = scb_mean( Y, level=level, method=method, param_method = param_method )$scb
      ### Check whether the SCB covers the true mean or respective the smoothed mean everywhere
      if( all(scb$lo<=mu.check) & all(scb$up>=mu.check) ){
        rate <- rate + 1
      }
    }
  }

  #### data.frame() returning the parameters and covering rate
  if( !is.null(param_scenario$SmoothWeights) ){
    retFrame <- list(covRate = rate/trials, covRateSmoothed = rate_smoothed/trials )
  }else{
    retFrame <- rate/trials
  }
  retFrame
}


#' Computing the covering rate of a given method for constructing simultaneous confidance
#' bands for the mean for a given error process.
#'
#' @param trials Integer amount of trials to obtain the covering rate.
#' @param N1 Integer sample size for sample 1.
#' @param N2 Integer sample size for sample 2.
#' @param level Numeric the targeted covering probability. Must be strictly between 0 and 1.
#' @param method String the name of the method used for quantile estimation. Possibilities are "tGKF", "GKF". Default value is "tGKF".
#' @param Mboots Numeric the number of bootstrap replicats used, if the quantile is estimated using a bootstrap method. Default value is 5000.
#' @param noise Function a function generating random realizations.
#' @param mean Vector a vector containing the values of the true mean of the random functions.
#' @return list with elements
#'  \itemize{
#'   \item cov.rate Numeric simulated covering rate
#'   \item trials Integer number of simulations to estimate the covering rate
#'   \item BootsReplicats Integer number of bootstrap replicates
#'   \item nSamp Integer number of
#' }
#' @export
covRate_2sample <- function(
  trials, N1 = 10, N2 = 10, level = 0.95, x = seq(0,1,length.out=100),
  mu1 = function(x) {sin(8*pi*x)*exp(-3*x)},
  mu2 = function(x) {sin(8*pi*x)*exp(-3*x)},
  noise1 = BernsteinSumNoise, noise2 = BernsteinSumNoise,  sigma1 = function(x){ ((0.6-x)^2+1)/6 },  sigma2 = function(x){ ((0.6-x)^2+1)/6 },
  sd_ObsNoise1 = 0, sd_ObsNoise2 = 0, method="tGKF", param_method = NULL,...
){

  mu1.sim = mu1(x)
  mu2.sim = mu2(x)

  rate  <- 0

  for(i in 1:trials){
    #### Generate data
    Y1 <- FunctionalDataSample( N=N1, x=x, mu=mu1, noise=noise1, sigma=sigma1,
                                sd_ObsNoise=sd_ObsNoise1,... ) ;
    Y2 <- FunctionalDataSample( N=N2, x=x, mu=mu2, noise=noise2, sigma=sigma2,
                                sd_ObsNoise=sd_ObsNoise2,... ) ;

    # #### PlotData
    # plot(NULL, xlim=c(0,1),ylim=c(min(y)-0.3, max(y)+0.3), xlab="Input", ylab="Response")
    # matlines(x,t(y), lty=1, col="lightblue")

    #### Confidence bands
    scb <- scb_meandiff( Y1, Y2, level = level, method=method, param_method=param_method )$scb

    ##### Check whether truth is covered
    if( all( 0 <= (scb$up )) & all( 0 >= (scb$lo )) ){
      rate <- rate + 1
    }
    #   matlines(xeval, t(rbind(mu.check, scb$hatmean, scb$scb)), lty=c(1,1,2,2), col=c("red","orange","black","black", "brown", "brown", "blue", "blue") )

  }

  #### data.frame() returning the parameters and covering rate
  data.frame( N1 = N1, N2 = N2, K = length(x), method = method, covRate = rate/trials )
}
