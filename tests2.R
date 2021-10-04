require(SIRF)
## test the generation of samples
M = FunctionalDataSample( N = 10,
                          x = seq( 0, 1, length.out = 100 ),
                          mu = function(x){ x^5 },
                          noise = SinCosSumNoise,
                          sigma = function(x){1},
                          sd_ObsNoise = 0 )

## test the covering rate function in 1D
covRate_simulation(
  trials = 2e3,
  scenario = "SpN",
  method = "tGKF",
  level = 0.95,
  N = 20,
  mu = function(x){0},
  noise = GaussDensitySumNoise,
  sigma = function(x) 1,
  sd_ObsNoise = 0,
  x = seq( 0, 1, length.out = 100 )
)

## test the covering rate function with smoothing in 1D
param_scenario <- list()
x                            = seq(0,1,length.out=100)
param_scenario$xeval         = seq(0,1,length.out=200)
param_scenario$SmoothWeights = locpol::locLinWeightsC( x=x, xeval=seq(0,1,length.out=200), bw=0.01, kernel=locpol::gaussK)$locWeig
covRate_simulation(
  trials   = 5000,
  scenario = "SpN",
  scenario = param_scenario,
  method   = "tGKF",
  level    = 0.95,
  N        = 10,
  mu       = function(x){ rep( 0, length(x) ) },
  noise    = GaussDensitySumNoise,
  sigma    = function(x) rep(1,length(x)), sd_ObsNoise=0, x=x
)

## test the covering rate function in 2D
N=10

covRate_simulation(
  trials   = 500,
  scenario = "SpN",
  method   = "tGKF", level=0.9, N=N, mu=function(x){ 3*outer( x, x, "*") }, noise=GaussKernelSum2DNoise, sigma=function(x) outer( x, x, FUN = function(s,t) (s+1)/(t^2+1) ), sd_ObsNoise=0, x=seq(0,1,length.out=100)
)

## test scale space
param_scenario        = list()
param_scenario$bw     = exp( seq( log(.02), log(0.1),length.out=50 ) )
kernel                = locpol::gaussK
param_scenario$kernel = kernel
param_scenario$xeval  = seq(0,1,length.out=100)

covRate_simulation(
  trials = 20, scenario = "Scale1d", param_scenario=param_scenario, method = "tGKF", level=0.20, N=10, mu=function(x){sin(x)}, noise=GaussDensitySumNoise, sigma=function(x) rep(1,length(x)), sd_ObsNoise=0, x=seq(0,1,length.out=100)
)


M = FunctionalDataSample( N=10, x=seq(0,1, length.out=50), mu=function(x){ 3*array(rep(array(outer( x, x, "*"), c(rep(length(x),2), 1)), N), c(rep(length(x),2), N))}, noise=toyNoise2Gauss2D, sigma=function(x) array(outer( x, x, FUN = function(s,t) (s+1)/(t^2+1) ), c(rep(length(x),2), 1)), sd_ObsNoise=0 )
