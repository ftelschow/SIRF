################################################################################################
##                                                                                          ####
##      Script providing the simulations for SCB paper                                      ####
##                                                                                          ####
################################################################################################
##
##
##
##
################################################################################################
#################################### Load required packages ####################################
################################################################################################
library(SCBfda)
library(fda)
rm(list=ls())

setwd("/media/sf_Linux/Research/Projects/2017_GKFinFDA")
################################################################################################
################################## one sample SCB simulations ##################################
################################################################################################
##
## Simulations for Section 8 of Telschow Schwartzman 2018
##
################################################################################################
rm(list=ls())
############ General parameters
x           = seq( 0, 1, length.out=100 )
level       = 0.95
trials      = 5e3
NVec        = c( 10, 20, 30, 50, 100, 150 )
sd_ObsNoise = 0
methodVec   = c("tGKF", "Bootstrapt", "MultiplierBootstrapt")
path_sim    = "/media/sf_Linux/Research/Projects/2017_GKFinFDA/Simulations/"
date        = "08_17_2018"

#################################### Simulations ####################################
################## Smooth Gaussian case
#### Model A
## Model specific parameters
mu_ModelA    = function(x){ sin(8*pi*x) * exp(-3*x) }
noise_ModelA = BernsteinSumNoise
sigma_ModelA = function(x){ ((1-x-0.4)^2+1)/6 }

## initialize the matrix for the covering rates
covRate_ModelAGauss <- matrix( NA, length(NVec), length(methodVec) )

## Simulation of the covering rates
for(k in 1:length(methodVec)){for( i in 1:length(NVec) ){
  covRate_ModelAGauss[i,k] = covRate_simulation( trials, scenario="SpN", method=methodVec[k], level=level, N=NVec[i], mu=mu_ModelA, noise=noise_ModelA, sigma=sigma_ModelA, sd_ObsNoise=0, x=x )[1,6]
}}

#### Model B
## Model specific parameters
mu_ModelB    = function(x){ sin(8*pi*x) * exp(-3*x) }
noise_ModelB = GaussDensitySumNoise
sigma_ModelB = function(x){ ((1-x-0.4)^2+1)/6 }

## initialize the matrix for the covering rates
covRate_ModelBGauss <- matrix( NA, length(NVec), length(methodVec) )

## Simulation of the covering rates
for(k in 1:length(methodVec)){for( i in 1:length(NVec) ){
  covRate_ModelBGauss[i,k] = covRate_simulation( trials, scenario="SpN", method=methodVec[k], level=level, N=NVec[i], mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=0, x=x )[1,6]
}}

#### Model C
## Model specific parameters
mu_ModelC    = function(x){ 3*outer( x, x, "*") }
noise_ModelC = GaussDensitySum2DNoise
sigma_ModelC = function(x){ outer( x, x, FUN = function(s,t) (s+1)/(t^2+1) ) }

## initialize the matrix for the covering rates
covRate_ModelCGauss <- matrix( NA, length(NVec), length(methodVec) )

## Simulation of the covering rates
for(k in 1:length(methodVec)){for( i in 1:length(NVec) ){
  covRate_ModelCGauss[i,k] = covRate_simulation( trials, scenario="SpN", method=methodVec[k], level=level, N=NVec[i], mu=mu_ModelC, noise=noise_ModelC, sigma=sigma_ModelC, sd_ObsNoise=0, x=seq(0, 1, length.out=50) )[1,6]
}}

## save the simulation
rm( k, i )
save.image( paste(path_sim, date,"_SimulationSmoothGauss.RData") )
rm( covRate_ModelAGauss, covRate_ModelBGauss, covRate_ModelCGauss )

################## Smooth Nongaussian case
#### Model A
## Model specific parameters
mu_ModelA    = function(x){ sin(8*pi*x) * exp(-3*x) }
noise_ModelA = function(N, x, sigma){ BernsteinSumNoise( N=N, x=x, sigma=sigma, randNumber=function(n){rt( n, df=3)/sqrt( 3/ (3-2) )} ) }
sigma_ModelA = function(x){ ((1-x-0.4)^2+1)/6 }

## initialize the matrix for the covering rates
covRate_ModelAnonGauss <- matrix( NA, length(NVec), length(methodVec) )

## Simulation of the covering rates
for(k in 1:length(methodVec)){for( i in 1:length(NVec) ){
  covRate_ModelAnonGauss[i,k] = covRate_simulation( trials, scenario="SpN", method=methodVec[k], level=level, N=NVec[i], mu=mu_ModelA, noise=noise_ModelA, sigma=sigma_ModelA, sd_ObsNoise=0, x=x )[1,6]
}}

#### Model B
## Model specific parameters
mu_ModelB    = function(x){ sin(8*pi*x) * exp(-3*x) }
noise_ModelB = function(N, x, sigma=function(x){ rep(1, length(x)) }, df){
                          randNumber=function(n){ (rchisq(n, df=df)-df)/sqrt(2*df)}
                          GaussDensitySumNoise( N=N, x=x, sigma=sigma, randNumber=randNumber ) }
sigma_ModelB = function(x){ ((1-x-0.4)^2+1)/6 }
dfVec        = c(1,2,7,10,30,60)

## initialize the matrix for the covering rates
covRate_ModelBnonGauss_tGKF <- covRate_ModelBnonGauss_boot <- matrix( NA, length(NVec), length(dfVec) )

## Simulation of the covering rates
for(k in 1:length(dfVec)){for( i in 1:length(NVec) ){
  covRate_ModelBnonGauss_tGKF[i,k] = covRate_simulation( trials, scenario="SpN", method="tGKF", level=level, N=NVec[i], mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=0, x=x, df=dfVec[k] )[1,6]
}}
for(k in 1:length(dfVec)){for( i in 1:length(NVec) ){
  covRate_ModelBnonGauss_boot[i,k] = covRate_simulation( trials, scenario="SpN", method="Bootstrapt", level=level, N=NVec[i], mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=0, x=x, df=dfVec[k] )[1,6]
}}

## save the simulation
rm(k, i)
save.image(paste(path_sim, date,"_SimulationSmoothNonGauss.RData"))
rm( covRate_ModelBnonGauss_tGKF, covRate_ModelBnonGauss_boot, covRate_ModelAnonGauss )

# ######################### Width mean and varaiance for different methods #########################
# # Compute the T statistic
# Tstat <- function(X, true.mu){
#   sqrt(dim(X)[2])*(rowMeans(X)-true.mu)/sqrt(matrixStats::rowVars(X))
# }
# # Simulation Parameter
# MonteN        = 1e5
# x             = seq(0,1, length.out=200)
# NSampVec      = length(NVec)
# Msim          = 1e3
# ################## Gaussian case
# # Parametsers
# noise    = GaussDensitySumNoise
#
# ############# Monte Carlo Simulation for the Correct Quantiles
# Is <- Sys.time()
# MonteTgauss <- matrix( NA, MonteN, NSampVec )
# for(k in 1:MonteN){
#   y = noise( N=max(NVec), x=x )
#   for(l in 1:NSampVec){
#     MonteTgauss[k,l] <- max(abs(Tstat( X=y[,1:NVec[l]], true.mu=rep(0, length(x)) )))
#   }
# }
# Ie <- Sys.time()
# Ie-Is
# rm(l,k)
#
# true.quantgauss <- apply(MonteTgauss, 2, function(vec) quantile(vec,level, type=8))
# print(true.quantgauss)
#
# hatquantilesgauss <- array( NA, dim=c( Msim, NSampVec, 7) )
#
# Is<-Sys.time()
# # Simulate the Distribution of the Quantiles
# for( m in 1:Msim ){
#   # Simulate the data
#   y = noise( N=max(NVec), x=x )
#   for( n in 1:NSampVec ){
#     # Generate data samples
#     data = y[,sample( 1:max(NVec), size = NVec[n] )]
#     # Estimate the LKCs
#     hatLKC = c(1, LKC_estim_direct( data ) )
#     # Compute the Quantiles of SCBs
#     hatquantilesgauss[m,n,1]  <- GKFquantileApprox( alpha = (1-level)/2, LKC=hatLKC, field="t", df=NVec[n] )
#     hatquantilesgauss[m,n,2]  <- GKFquantileApprox( alpha = (1-level)/2, LKC=hatLKC, field="Gauss" )
#     degras <- SCBmeanfd::scb.mean( x, t(data), bandwidth=diff(x)[1]/3, level=level, scbtype="normal", nboot=5e3 )
#     hatquantilesgauss[m,n,3]  <- degras$qnorm
#     hatquantilesgauss[m,n,4]  <- NonParametricBootstrap( A=data, params=list(Mboots=5e3, alpha=1-level, method="regular") )$q
#     hatquantilesgauss[m,n,5]  <- NonParametricBootstrap( A=data, params=list(Mboots=5e3, alpha=1-level, method="t") )$q
#     hatquantilesgauss[m,n,6]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="regular") )$q
#     hatquantilesgauss[m,n,7]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="t") )$q
#   }
# }
# Ie <- Sys.time()
# SimulationTimeQuant = Ie-Is
# rm(n,m, Ie, Is, y, data)
#
# ################## non Gaussian case
# # Parametsers
# noise    = function(N, x, sigma=function(x){ rep(1, length(x)) }, df=7){
#   randNumber=function(n){ (rchisq(n, df=df)-df)/sqrt(2*df)}
#   GaussDensitySumNoise( N=N, x=x, sigma=sigma, randNumber=randNumber ) }
#
# ############# Monte Carlo Simulation for the Correct Quantiles
# Is <- Sys.time()
# MonteTnongauss <- matrix( NA, MonteN, NSampVec )
# for(k in 1:MonteN){
#   y = noise( N=max(NVec), x=x )
#   for(l in 1:NSampVec){
#     MonteTnongauss[k,l] <- max(abs(Tstat( X=y[,1:NVec[l]], true.mu=rep(0, length(x)) )))
#   }
# }
# Ie <- Sys.time()
# Ie-Is
# rm(l,k)
#
# true.quantnongauss <- apply(MonteTnongauss, 2, function(vec) quantile(vec,level, type=8))
# print(true.quantnongauss)
#
# hatquantilesnongauss <- array( NA, dim=c( Msim, NSampVec, 7) )
#
# Is<-Sys.time()
# # Simulate the Distribution of the Quantiles
# for( m in 1:Msim ){
#   # Simulate the data
#   y = noise( N=max(NVec), x=x )
#   for( n in 1:NSampVec ){
#     # Generate data samples
#     data = y[,sample( 1:max(NVec), size = NVec[n] )]
#     # Estimate the LKCs
#     hatLKC = c(1, LKC_estim_direct( data ) )
#     # Compute the Quantiles of SCBs
#     hatquantilesnongauss[m,n,1]  <- GKFquantileApprox( alpha = (1-level)/2, LKC=hatLKC, field="t", df=NVec[n] )
#     hatquantilesnongauss[m,n,2]  <- GKFquantileApprox( alpha = (1-level)/2, LKC=hatLKC, field="Gauss" )
#     degras <- SCBmeanfd::scb.mean( x, t(data), bandwidth=diff(x)[1]/3, level=level, scbtype="normal", nboot=5e3 )
#     hatquantilesnongauss[m,n,3]  <- degras$qnorm
#     hatquantilesnongauss[m,n,4]  <- NonParametricBootstrap( A=data, params=list(Mboots=5e3, alpha=1-level, method="regular") )$q
#     hatquantilesnongauss[m,n,5]  <- NonParametricBootstrap( A=data, params=list(Mboots=5e3, alpha=1-level, method="t") )$q
#     hatquantilesnongauss[m,n,6]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="regular") )$q
#     hatquantilesnongauss[m,n,7]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="t") )$q
#   }
# }
# Ie <- Sys.time()
# SimulationTimeQuantnon = Ie-Is
# rm(n,m, Ie, Is, y, data)
#
# ## Save Workspace
# save.image(paste(path_sim, date,"_SimulationQuantiles.RData",sep="") )
# rm(hatquantilesnongauss, hatquantilesgauss, MonteTnongauss, true.quantnongauss, MonteTgauss, true.quantgauss )

################## Observation and bandwidth dependence
#### Parameters
obsSDvec = c(0.05, 0.1, 0.2)
BWvec    = c(0.02, 0.05, 0.1)
#### Model A
## Model specific parameters
mu_ModelA    = function(x){ sin(8*pi*x) * exp(-3*x) }
noise_ModelA = BernsteinSumNoise
sigma_ModelA = function(x){ ((1-x-0.4)^2+1)/6 }

## initialize the matrix for the covering rates
covRate_ModelANonSmooth <- array( NA, dim = c(length(obsSDvec), length(NVec), length(BWvec)) )

## Simulation of the covering rates
for( k in 1:length(BWvec) ){
  # precompute the smoothing matrix
  param_scenario <- list()
  param_scenario$xeval         = seq(0,1,length.out=300)
  param_scenario$SmoothWeights = locpol::locLinWeightsC( x=x, xeval=param_scenario$xeval , bw=BWvec[k], kernel=locpol::gaussK)$locWeig
  for( j in 1:length(NVec) ){for(i in 1:length(obsSDvec)){
    param_scenario$SmoothWeights
    covRate_ModelANonSmooth[i,j,k] = covRate_simulation( trials, scenario="SpN", param_scenario = param_scenario, method="tGKF", level=level, N=NVec[i], mu=mu_ModelA, noise=noise_ModelA, sigma=sigma_ModelA, sd_ObsNoise=obsSDvec[i], x=x )[1,7]
  }}
}

#### Model B
## Model specific parameters
mu_ModelB    = function(x){ sin(8*pi*x) * exp(-3*x) }
noise_ModelB = GaussDensitySumNoise
sigma_ModelB = function(x){ ((1-x-0.4)^2+1)/6 }

## initialize the matrix for the covering rates
covRate_ModelBNonSmooth <- array( NA, dim = c(length(obsSDvec), length(NVec), length(BWvec)) )

## Simulation of the covering rates
for( k in 1:length(BWvec) ){
  # precompute the smoothing matrix
  param_scenario <- list()
  param_scenario$xeval         = seq(0,1,length.out=300)
  param_scenario$SmoothWeights = locpol::locLinWeightsC( x=x, xeval=param_scenario$xeval, bw=BWvec[k], kernel=locpol::gaussK)$locWeig
  for( j in 1:length(NVec) ){for(i in 1:length(obsSDvec)){
    param_scenario$SmoothWeights
    covRate_ModelBNonSmooth[i,j,k] = covRate_simulation( trials, scenario="SpN", param_scenario = param_scenario, method="tGKF", level=level, N=NVec[i], mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=obsSDvec[i], x=x )[1,7]
  }}
}
#### Scale space
## initialize the matrix for the covering rates
covRate_ModelBScale <- array( NA, dim = c(length(NVec), length(obsSDvec)) )

# precompute the smoothing matrix
param_scenario        = list()
param_scenario$bw     = exp( seq( log(.02), log(0.1),length.out=100 ) )
kernel                = locpol::gaussK
param_scenario$kernel = kernel
param_scenario$xeval  = seq(0,1,length.out=100)

## Simulation of the covering rates
for( j in 1:length(NVec) ){for(i in 1:length(obsSDvec)){
  param_scenario$SmoothWeights
  covRate_ModelBScale[j,i] = covRate_simulation( trials, scenario="Scale1d", param_scenario = param_scenario, method="tGKF", level=level, N=NVec[i], mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=obsSDvec[i], x=x )[1,6]
}}

## save the simulation
rm(k, i, j)
save.image( paste(path_sim, date,"_SimulationNonSmoothGauss.RData") )
rm(covRate_ModelANonSmooth, covRate_ModelBNonSmooth, covRate_ModelBScale)

# ############ Plots of the results
# clab  = 1.5
# caxis = 1.5
# cmain = 2
#
# ###### Plot the results of the smooth gaussian one sample case
# load("/media/sf_Linux/Research/Projects/2017_GKFinFDA/Workspaces/ArticleSimulations/20180807_SimulationSmoothGauss.RData")
# pdfname <- "Pics/ResultsSimulationSmoothGaussian.pdf"
# pdf(pdfname,  title=pdfname, width=1.1*10, height=1.1*7.5)
# par(mfrow=c(3,3))
# par(mar=c(2.1,5.1,3.1,2.1) )
# #### Plot Sample paths of Data
# set.seed(1)
# # Model A plot
# y = FunctionalDataSample( N=15, x=x, mu=mu_ModelA, noise=noise_ModelA, sigma=sigma_ModelA )
# plot( NULL, xlim=c(0,1), ylim=range(y), xlab="", ylab="", main="Model A", cex.main=cmain, cex.axis=caxis )
# matlines(x, y)
# lines(x, mu_ModelA(x), lwd=2)
#
# # Model B plot
# y = FunctionalDataSample( N=15, x=x, mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB )
# plot( NULL, xlim=c(0,1), ylim=range(y), xlab="", ylab="", main="Model B", cex.main=cmain, cex.axis=caxis )
# matlines(x, y)
# lines(x, mu_ModelB(x), lwd=2)
#
# # Model C plot
# y = FunctionalDataSample( N=15, x=seq(0,1,length.out=50), mu=mu_ModelC, noise=noise_ModelC, sigma=sigma_ModelC )
# fields::image.plot(y[,,9], main="Model C", cex.axis=caxis, cex.main=cmain)
#
# #### Plot Sample paths of Noise
# # Noise Model A plot
# y = noise_ModelA( N=15, x=x, sigma=function(x){rep(1, length(x))} )
# plot(NULL, xlim=c(0,1), ylim=range(y), xlab="", ylab="", main="Noise Model A", cex.main=cmain, cex.axis=caxis)
# matlines(x, y)
#
# # Model B plot
# y = noise_ModelB( N=15, x=x, sigma=function(x){rep(1, length(x))} )
# plot(NULL, xlim=c(0,1), ylim=range(y), xlab="", ylab="", main="Noise Model B", cex.main=cmain, cex.axis=caxis)
# matlines(x, y)
#
# # Model C plot
# y = noise_ModelC( N=15, x=x )
# fields::image.plot(y[,,9], main="Noise Model C", cex.axis=caxis, cex.main=cmain)
#
#
# par(mar=c(5.1,5.1,2.1,2.1) )
# ### Plot Results
# ## constants
# colVec <- c("red", "blue", "cyan3")
# axisPoints = 10*(1:length(NVec))
#
# ## Model A results
# # setup the plot
# plot( NULL, xlim = range(axisPoints), ylim = c(86,100), ylab="Covering Rate [%]", xlab="Sample Size [N]", xaxt='n', cex.lab=clab, cex.axis=caxis )
# # set axis labels
# axis(1, at=axisPoints, labels=NVec, cex.axis=caxis)
# # nominal level line and confidence regions
# abline(h=level*100, col="black", lty=1)
# abline( h = 100*(level+c(-2,2)*sqrt(level*(1-level)/trials)), col="black", lty=2)
# # Plot the data
# matlines( axisPoints, 100*covRate_ModelAGauss, col=colVec, lty=rep(1, length(colVec)) )
# legend("bottomright", legend=c("NaiveBootstrap", "tGKF", "MultiplierBootstrap"), col=colVec[c(2,1,3)], pch=c(2,1,3), lty=rep(1,3), cex=1.5)
#
# ## Model B results
# # setup the plot
# plot( NULL, xlim = range(axisPoints), ylim = c(86,100), ylab="Covering Rate [%]", xlab="Sample Size [N]", xaxt='n', cex.lab=clab, cex.axis=caxis )
# # set axis labels
# axis(1, at=axisPoints, labels=NVec, cex.axis=caxis)
# # nominal level line and confidence regions
# abline(h=level*100, col="black", lty=1)
# abline( h = 100*(level+c(-2,2)*sqrt(level*(1-level)/trials)), col="black", lty=2)
# # Plot the data
# matlines( axisPoints, 100*covRate_ModelBGauss, col=colVec, lty=rep(1, length(colVec)) )
# legend("bottomright", legend=c("NaiveBootstrap", "tGKF", "MultiplierBootstrap"), col=colVec[c(2,1,3)], pch=c(2,1,3), lty=rep(1,3), cex=1.5)
#
# ## Model C results
# # setup the plot
# plot( NULL, xlim = range(axisPoints), ylim = c(86,100), ylab="Covering Rate [%]", xlab="Sample Size [N]", xaxt='n', cex.lab=clab, cex.axis=caxis )
# # set axis labels
# axis(1, at=axisPoints, labels=NVec, cex.axis=caxis)
# # nominal level line and confidence regions
# abline(h=level*100, col="black", lty=1)
# abline( h = 100*(level+c(-2,2)*sqrt(level*(1-level)/trials)), col="black", lty=2)
# # Plot the data
# matlines( axisPoints, 100*covRate_ModelCGauss, col=colVec, lty=rep(1, length(colVec)) )
# legend("bottomright", legend=c("NaiveBootstrap", "tGKF", "MultiplierBootstrap"), col=colVec[c(2,1,3)], pch=c(2,1,3), lty=rep(1,3), cex=1.5)
#
# dev.off()
#
# ###### Plot the results of the smooth Non-gaussian one sample case
#
# ################################################################################################
# ################################## two sample SCB simulations ##################################
# ################################################################################################
# ##
# ## Simulations for Section ?? of Telschow Schwartzman 2018
# ##
# ################################################################################################
# rm(list=ls())
# ############ General parameters
# x           = seq(0,1,length.out=100)
# level       = 0.95
# trials      = 5e3
# N           = cbind( c(10, 10), c(20, 20), c(30, 30),  c(50, 50),  c(100, 100), c(150, 150), c(20, 40), c(80, 160), c(150, 300), c(20, 80), c(50, 200), c(100, 400)  )
# sd_ObsNoise = 0
# I           = dim(N)[2]
#
# path_sim    = "/media/sf_Linux/Research/Projects/2017_GKFinFDA/Simulations/"
# covRate_SameCorSameVar <- covRate_DiffCorSameVar <- covRate_SameCorDiffVar <- covRate_DiffCorDiffVar <- NULL
# ############ Simulations
# ## same correlation, same variance
# noise1 = BernsteinSumNoise
# noise2 = BernsteinSumNoise
# sigma1 = function(x){ ((0.6-x)^2+1)/6 }
# sigma2 = function(x){ ((0.6-x)^2+1)/6 }
#
# for( i in 1:I ){
#   covRate_SameCorSameVar <- rbind( covRate_SameCorSameVar, covRate_2sample(
#     trials, N1 = N[1,i], N2 = N[2,i], level = level, x = x,
#     noise1 = noise1, noise2 = noise2,  sigma1 = sigma1,  sigma2 = sigma2
#   ))
#   print(covRate_SameCorSameVar)
# }
#
#
# ## different correlation, same variance
# noise1 = BernsteinSumNoise
# noise2 = GaussDensitySumNoise
# sigma1 = function(x){ ((0.6-x)^2+1)/6 }
# sigma2 = function(x){ ((0.6-x)^2+1)/6 }
#
# for( i in 1:I ){
#   covRate_DiffCorSameVar <- rbind( covRate_DiffCorSameVar, covRate_2sample(
#     trials, N1 = N[1,i], N2 = N[2,i], level = level, x = x,
#     noise1 = noise1, noise2 = noise2,  sigma1 = sigma1,  sigma2 = sigma2
#   ))
#   print(covRate_DiffCorSameVar)
# }
#
# ## same covariance, different variance
# noise1 = BernsteinSumNoise
# noise2 = BernsteinSumNoise
# sigma1 = function(x){ ((0.6-x)^2+1)/6 }
# sigma2 = function(x){ rep(0.2, length(x)) }
#
# for( i in 1:I ){
#   covRate_SameCorDiffVar <- rbind( covRate_SameCorDiffVar, covRate_2sample(
#     trials, N1 = N[1,i], N2 = N[2,i], level = level, x = x,
#     noise1 = noise1, noise2 = noise2,  sigma1 = sigma1,  sigma2 = sigma2
#   ))
#   print(covRate_SameCorDiffVar)
# }
#
# ## different correlation, different variance
# noise1 = BernsteinSumNoise
# noise2 = GaussDensitySumNoise
# sigma1 = function(x){ ((0.6-x)^2+1)/6 }
# sigma2 = function(x){ rep(0.2, length(x)) }
#
# for( i in 1:I ){
#   covRate_DiffCorDiffVar <- rbind( covRate_DiffCorDiffVar, covRate_2sample(
#     trials, N1 = N[1,i], N2 = N[2,i], level = level, x = x,
#     noise1 = noise1, noise2 = noise2,  sigma1 = sigma1,  sigma2 = sigma2
#   ))
#   print(covRate_DiffCorDiffVar)
# }
#
# rm(i, I, noise1, noise2, sigma1, sigma2 )
# # save.image(paste(path_sim, "Simulation_GKFinFDA_SCBdiffmean.RData"))
#
# ############ Plot Simulations Results
# clab  = 1.5
# caxis = 1.5
# cmain = 2
#
# ###### Plot the results of the smooth gaussian one sample case
# load(paste(path_sim, "Simulation_GKFinFDA_SCBdiffmean.RData"))
# pdfname <- "Pics/ResultsSimulationSCBdiffmean.pdf"
# pdf(pdfname,  title=pdfname, width=1.45*10, height=1.1*3.5)
# par(mfrow=c(1,3))
# par(mar=c(5.1,5.1,3.1,2.1) )
#   Nvec   = (1:6)*20#covRate_SameCorSameVar[1:6,1]
#   colvec = c("cyan3", "red", "blue", "darkgreen")
#   pchvec = 1:4
#
#   ## Plot balanced sample size
#   plot( NULL, xlim=range(Nvec), main="c=1", ylim=100*c(.925,0.96), ylab="Covering Rate [%]", xlab="Sample Size [N]", xaxt='n', cex.lab=clab, cex.axis=caxis, cex.main=cmain)
#   axis(1, at=Nvec, labels=c("10/10", "20/20", "30/30", "50/50", "100/100", "150/150"), cex.axis=caxis)
#   axis(1, at=Nvec[6], labels="150/150", cex.axis=caxis)
#   abline( h = 100*0.95, col="black" )
#   abline( h = 100*0.95+100*c(1,-1)*2*sqrt(0.95*0.05/trials), col="black", lty=2)
#
#   lines(  Nvec, 100*covRate_SameCorSameVar[1:6,5], col=colvec[1] )
#   points( Nvec, 100*covRate_SameCorSameVar[1:6,5], col=colvec[1], pch=pchvec[1] )
#   lines(  Nvec, 100*covRate_DiffCorSameVar[1:6,5], col=colvec[2])
#   points( Nvec, 100*covRate_DiffCorSameVar[1:6,5], col=colvec[2], pch=pchvec[2] )
#   lines(  Nvec, 100*covRate_SameCorDiffVar[1:6,5], col=colvec[3])
#   points( Nvec, 100*covRate_SameCorDiffVar[1:6,5], col=colvec[3], pch=pchvec[3] )
#   lines(  Nvec, 100*covRate_DiffCorDiffVar[1:6,5], col=colvec[4])
#   points( Nvec, 100*covRate_DiffCorDiffVar[1:6,5], col=colvec[4], pch=pchvec[4] )
#   legend("bottomright", legend=c("Same Corr, Same Var", "Diff Corr, Same Var", "Same Corr, Diff Var", "Diff Corr, Diff Var"), col=colvec, pch=pchvec, lty=rep(1,length(colvec)), cex=1.5)
#
#   ## plot Factor 2 between sample size
#   Nvec   = 20*(1:3) #covRate_SameCorSameVar[7:9,1]
#
#   plot( NULL, xlim=range(Nvec), main="c=2", ylim=100*c(.925,0.96), ylab="Covering Rate [%]", xlab="Sample Size [N]", xaxt='n', cex.lab=clab, cex.axis=caxis, cex.main=cmain)
#   axis(1, at=Nvec, labels=c("20/40", "80/160", "150/300"), cex.axis=caxis)
# #  axis(1, at=Nvec[6], labels="150/150", cex.axis=caxis)
#
#   abline( h = 95, col="black" )
#   abline( h = 95+100*c(1,-1)*2*sqrt(0.95*0.05/trials), col="black", lty=2)
#
#   lines(  Nvec, 100*covRate_SameCorSameVar[7:9,5], col=colvec[1] )
#   points( Nvec, 100*covRate_SameCorSameVar[7:9,5], col=colvec[1], pch=pchvec[1] )
#   lines(  Nvec, 100*covRate_DiffCorSameVar[7:9,5], col=colvec[2])
#   points( Nvec, 100*covRate_DiffCorSameVar[7:9,5], col=colvec[2], pch=pchvec[2] )
#   lines(  Nvec, 100*covRate_SameCorDiffVar[7:9,5], col=colvec[3])
#   points( Nvec, 100*covRate_SameCorDiffVar[7:9,5], col=colvec[3], pch=pchvec[3] )
#   lines(  Nvec, 100*covRate_DiffCorDiffVar[7:9,5], col=colvec[4])
#   points( Nvec, 100*covRate_DiffCorDiffVar[7:9,5], col=colvec[4], pch=pchvec[4] )
#
#
#   ## plot factor 4 between sample size
#   Nvec   = 1:3*20#covRate_SameCorSameVar[10:12,1]
#
#   plot( NULL, xlim=range(Nvec), main="c=4", ylim=100*c(.925,0.96), ylab="Covering Rate [%]", xlab="Sample Size [N]", xaxt='n', cex.lab=clab, cex.axis=caxis, cex.main=cmain)
#   axis(1, at=Nvec, labels=c("20/80", "50/200", "100/400"), cex.axis=caxis)
#   abline( h = 95, col="black" )
#   abline( h = 95+100*c(1,-1)*2*sqrt(0.95*0.05/trials), col="black", lty=2)
#
#   lines(  Nvec, 100*covRate_SameCorSameVar[10:12,5], col=colvec[1] )
#   points( Nvec, 100*covRate_SameCorSameVar[10:12,5], col=colvec[1], pch=pchvec[1] )
#   lines(  Nvec, 100*covRate_DiffCorSameVar[10:12,5], col=colvec[2])
#   points( Nvec, 100*covRate_DiffCorSameVar[10:12,5], col=colvec[2], pch=pchvec[2] )
#   lines(  Nvec, 100*covRate_SameCorDiffVar[10:12,5], col=colvec[3])
#   points( Nvec, 100*covRate_SameCorDiffVar[10:12,5], col=colvec[3], pch=pchvec[3] )
#   lines(  Nvec, 100*covRate_DiffCorDiffVar[10:12,5], col=colvec[4])
#   points( Nvec, 100*covRate_DiffCorDiffVar[10:12,5], col=colvec[4], pch=pchvec[4] )
# dev.off()
