################################################################################################
##                                                                                          ####
##      Script providing the simulations for FTAS_2018_GKFinFDA.tex                         ####
##                                                                                          ####
################################################################################################
#### Prepare work space
rm( list = ls() )
# Set working directory
WD  <- "/home/drtea/Research/Projects/2017_GKFinFDA"
setwd( WD )

path_sim  = paste( WD, "/Simulations/", sep = "" )
path_pics = paste( WD, "/Pics/", sep = "" )

#### Graphical parameters
clab  = 1.5
caxis = 1.5
cmain = 2

################################################################################################
#################################### Load required packages ####################################
################################################################################################
library( xtable )
library( data.table )
library( rmutil )
library( SCBfda )
library( fda )
library( ffscb )


################################################################################################
################################## one sample SCB simulations ##################################
################################################################################################
##
## Simulations for Section 8 of Telschow Schwartzman 2018
##
################################################################################################
################################################################################################
##                                                                                            ##
##                          Simulations and Plots Smooth Gaussian Case                        ##
##                                                                                            ##
################################################################################################
############ General parameters
x           = seq( 0, 1, length.out = 100 )
level       = 0.95
Msim        = 5e3
Mboots      = 5e3
NVec        = c( 10, 20, 30, 50, 100, 150 )
sd_ObsNoise = 0
methodVec   = c( "tGKF",
                 "NonParametricBootstrap",
                 "MultiplierBootstrap",
                 "MultiplierBootstrap",
                 "ffscb" )
param_method            = list()
param_method$tGKF       = list( L0 = 1, LKC_estim = LKC_estim_direct )
param_method$Bootstrapt = list( Mboots = Mboots )
param_method$gMultiplierBootstrapt = list( Mboots = Mboots, weights = "gauss" )
param_method$rMultiplierBootstrapt = list( Mboots = Mboots, weights = "rademacher" )
param_method$ffscb = list(  )

methodVec2D    = methodVec[ 1:4 ]
param_method2D = param_method[ 1:4 ]

date        = "y2021_m01_d06"

############ Simulations
#### Model A
# Model specific parameters
mu_ModelA    = function( x ){ sin( 8 * pi * x ) * exp( -3 * x ) }
noise_ModelA = BernsteinSumNoise
sigma_ModelA = function(x){ ( ( 1 - x - 0.4 )^2 + 1 ) / 6 }

# Simulation of the covering rates
Is <- Sys.time()
covRate_ModelAGauss = t( covRate_simulation( Msim,
                                             scenario = "SpN",
                                             method   = methodVec,
                                             param_method = param_method,
                                             level = level,
                                             N     = NVec,
                                             mu    = mu_ModelA,
                                             noise = noise_ModelA,
                                             sigma = sigma_ModelA,
                                             sd_ObsNoise = 0, x = x ) )
Ie <- Sys.time()
# Show simulation time
simulationTimeModelAsmoothGauss <- Ie - Is
print( Ie - Is )

# Timing model A smooth Gauss
# Simulation of the covering rates
simTimesModelAsmoothGauss <- 1:length( methodVec )
names( simTimesModelAsmoothGauss ) <- methodVec

for( k in 1:length( methodVec ) ){
  Is <- proc.time()
  covRate_simulation( 2e1,
                      scenario = "SpN",
                      method   = methodVec[k],
                      param_method = param_method[k],
                      level = level,
                      N     = NVec,
                      mu    = mu_ModelA,
                      noise = noise_ModelA,
                      sigma = sigma_ModelA,
                      sd_ObsNoise = 0,
                      x     = x )
  Ie <- proc.time()
  # Show simulation time
  tmp = Ie - Is
  simTimesModelAsmoothGauss[k] <- tmp[3]
  print( Ie - Is )
}


#### Model B
# Model specific parameters
mu_ModelB    = function( x ){ sin( 8 * pi * x ) * exp( -3 * x ) }
noise_ModelB = GaussDensitySumNoise
sigma_ModelB = function( x ){ ( ( 1 - x - 0.4 )^2 + 1 ) / 6 }

# Simulation of the covering rates
Is <- Sys.time()
covRate_ModelBGauss = t( covRate_simulation( Msim,
                                             scenario = "SpN",
                                             method   = methodVec,
                                             param_method = param_method,
                                             level = level,
                                             N     = NVec,
                                             mu    = mu_ModelB,
                                             noise = noise_ModelB,
                                             sigma = sigma_ModelB,
                                             sd_ObsNoise = 0,
                                             x     = x ) )
Ie <- Sys.time()
# Show simulation time
simulationTimeModelBsmoothGauss <- Ie-Is
print(Ie-Is)

#### Model C
# Model specific parameters
mu_ModelC    = function( x ){ 3 * outer( x, x, "*" ) }
noise_ModelC = GaussDensitySum2DNoise
sigma_ModelC = function( x ){
                    outer( x, x, FUN = function( s, t ) ( s + 1 ) / ( t^2 + 1 ) )
                }

# Simulation of the covering rates
Is <- proc.time()
x2d = seq( 0, 1, length.out = 50)
covRate_ModelCGauss = t( covRate_simulation( Msim,
                                             scenario = "SpN",
                                             method   = methodVec2D,
                                             param_method = param_method2D,
                                             level = level,
                                             N     = NVec,
                                             mu    = mu_ModelC,
                                             noise = noise_ModelC,
                                             sigma = sigma_ModelC,
                                             sd_ObsNoise = 0,
                                             x     = x2d ) )
Ie <- proc.time()
# Show simulation time
simulationTimeModelCsmoothGauss <- Ie - Is
print( Ie - Is )

# Timing model C smooth Gauss
# Simulation of the covering rates
  simTimesModelCsmoothGauss <- 1:length( methodVec2D )
  names( simTimesModelCsmoothGauss ) <- methodVec2D

  for( k in 1:length( methodVec2D ) ){
    Is <- Sys.time()
    covRate_simulation( 2e1,
                        scenario = "SpN",
                        method   = methodVec2D[k],
                        param_method = param_method2D[k],
                        level = level,
                        N     = NVec,
                        mu    = mu_ModelC,
                        noise = noise_ModelC,
                        sigma = sigma_ModelC,
                        sd_ObsNoise = 0,
                        x     = x )
    Ie <- Sys.time()
    # Show simulation time
    simTimesModelCsmoothGauss[ k ] <- Ie - Is
    print( Ie - Is )
  }

  #### save the simulation
  rm( k, i, Ie, Is )
  save.image( paste(path_sim, date,"_SimulationSmoothGauss.RData") )

############ Plot the Results and error processes
date = "y2021_m01_d06"
load( paste( path_sim, date, "_SimulationSmoothGauss.RData" ) )
pdfname <- paste( path_pics, "ResultsSimulationSmoothGaussian.pdf", sep = "" )
pdf( pdfname,  title = pdfname, width = 1.1 * 10, height = 1.1 * 7.5 )
par( mfrow = c( 3, 3 ) )
  #### Plot Sample paths of Signal plus noise model
  set.seed(1)
  # Model A sample paths
  x = seq(0,1, length.out=100)
  y = FunctionalDataSample( N=15, x=x, mu=mu_ModelA, noise=noise_ModelA, sigma=sigma_ModelA, sd_ObsNoise=0 )
      par(mar=c(2.1,5.1,3.1,2.1) )
  plot(NULL, xlim=c(0,1), ylim=range(y), xlab="", ylab="", main="Model A", cex.main=cmain, cex.axis=caxis)
  matlines(x, y)
  lines(x, mu_ModelA(x), lwd=2)

  # Model B sample paths plot
  x=seq(0,1, length.out=100)
  y = FunctionalDataSample( N=15, x=x, mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=0 )
  par(mar=c(2.1,5.1,3.1,2.1) )
  plot(NULL, xlim=c(0,1), ylim=range(y), xlab="", ylab="", main="Model B", cex.main=cmain, cex.axis=caxis)
  matlines(x, y)
  lines(x, mu_ModelB(x), lwd=2)

  # Model C sample paths plot
  x = seq(0,1,length.out=50)
  y = FunctionalDataSample( N=15, x=x, mu=mu_ModelC, noise=noise_ModelC, sigma=sigma_ModelC, sd_ObsNoise=0 )
  par(mar=c(2.1,5.1,3.1,5) )
  fields::image.plot(y[,,9], main="Model C", cex.axis=caxis, cex.main=cmain)

  #### Plot sample paths of noise
  # Noise Model A plot
  x = seq( 0, 1, length.out = 100 )
  y = noise_ModelA( N = 15, x = x )
  par( mar = c( 2.1, 5.1, 3.1, 2.1 ) )
  plot( NULL, xlim = c( 0, 1 ), ylim = range( y ),
        xlab = "", ylab = "",
        main = "Noise Model A",
        cex.main = cmain, cex.axis = caxis )
  matlines(x, y)

  # Noise Model B plot
  x = seq( 0, 1, length.out = 100 )
  y = noise_ModelB( N = 15, x = x )
  par( mar = c( 2.1, 5.1, 3.1, 2.1 ) )
  plot( NULL, xlim = c( 0, 1 ), ylim = range( y ),
        xlab = "", ylab = "",
        main = "Noise Model B",
        cex.main = cmain, cex.axis = caxis )
  matlines(x, y)

  # Noise Model C plot
  x = seq( 0, 1, length.out = 50 )
  y = noise_ModelC( N = 10, x = x )
  par( mar = c( 2.1, 5.1, 3.1, 5 ) )
  fields::image.plot( y[ , , 9 ],
                      main = "Noise Model C",
                      cex.axis = caxis,
                      cex.main = cmain )
  # Clear Workspace
  rm(x, y)

  #### Plot Simulation Results
  par( mar = c( 5.1, 5.1, 2.1, 2.1 ) )
  colVec <- c("red", "darkmagenta", "cyan3", "blue", "orange")
  for(i in 1:3){
    if(i==1){
      cov.data = covRate_ModelAGauss
      model    = "A"
      Nmethods = length(colnames(cov.data))
    }else if(i==2){
      cov.data = covRate_ModelBGauss
      model    = "B"
    }else if(i==3){
      cov.data = covRate_ModelCGauss
      model    = "C"
    }
    axisPoints = 10*(1:length(NVec))
    plot( NULL, xlim = range(axisPoints), ylim = c(84,100), ylab="Covering Rate [%]", xlab="Number of Samples [N]", xaxt='n', cex.lab=clab, cex.axis=caxis )
    axis(1, at=axisPoints, labels=NVec, cex.axis=caxis)

    ### nominal level
    abline(h=level*100, col="black", lty=1)
    ### 2 times standard error of nominal level
    abline( h = 100*(level+c(-2,2)*sqrt(level*(1-level)/Msim)), col="black", lty=2)

    ### Plot the data
    Nlength = length(NVec)
    for(k in 1:Nmethods){
      lines(  axisPoints, 100 * cov.data[ ,k ], col = colVec[k], pch = k  )
      points( axisPoints, 100 * cov.data[ ,k ], col = colVec[k], pch = k  )
    }
    if(i==1){
      legend( "bottomright",
              legend = c( "Boots-t",
                          "tGKF",
                          "rMult-t",
                          "gMult-t",
                          "ffscb" ),
              col = colVec[ c( 2, 1, 4, 3, 5) ],
              pch = c( 2, 1, 4, 3, 5),
              lty = rep( 1, Nmethods ), cex = 1.2 )
    }
    rm(k, axisPoints, Nlength, model, cov.data)
  }
dev.off()

#### Clear Workspace
rm( covRate_ModelAGauss, covRate_ModelBGauss, covRate_ModelCGauss, i, colVec )
rm( x, level, Msim, NVec, sd_ObsNoise, methodVec, date )
rm( mu_ModelA, mu_ModelB, mu_ModelC, noise_ModelA,noise_ModelB, noise_ModelC, sigma_ModelA, sigma_ModelB, sigma_ModelC )
rm( simulationTimeModelAsmoothGauss, simulationTimeModelBsmoothGauss, simulationTimeModelCsmoothGauss )

################################################################################################
##                                                                                            ##
##                      Simulations and Plots Smooth Non-Gaussian Case                        ##
##                                                                                            ##
################################################################################################
############ General parameters
x           = seq( 0, 1, length.out=100 )
level       = 0.95
Msim        = 5e3
NVec        = c( 10, 20, 30, 50, 100, 150 )
sd_ObsNoise = 0
methodVec   = c( "tGKF",
                 "NonParametricBootstrap",
                 "MultiplierBootstrap",
                 "MultiplierBootstrap",
                 "ffscb" )
param_method = list()
param_method$tGKF       = list( L0=1, LKC_estim = LKC_estim_direct )
param_method$Bootstrapt = list( Mboots=5e3 )
param_method$gMultiplierBootstrapt = list( Mboots=5e3, weights="gauss")
param_method$rMultiplierBootstrapt = list( Mboots=5e3, weights="rademacher")
param_method$ffscb = list( )
date        = "11_27_2020"

############ Simulations
#### Model A
# Model specific parameters
mu_ModelA    = function(x){ sin( 8 * pi * x) * exp( -3 * x ) }
noise_ModelA = function(N, x, sigma){
            BernsteinSumNoise( N = N,
            x = x,
            sigma = sigma,
            randNumber = function( n ){ rt( n, df = 3 ) / sqrt( 3 / ( 3 - 2 ) )
              } ) }
sigma_ModelA = function(x){ ((1-x-0.4)^2+1)/6 }

# Simulation of the covering rates
Is <- Sys.time()
covRate_ModelAnonGauss = t(covRate_simulation( Msim, scenario="SpN", method=methodVec, param_method=param_method, level=level, N=NVec, mu=mu_ModelA, noise=noise_ModelA, sigma=sigma_ModelA, sd_ObsNoise=0, x=x ))
Ie <- Sys.time()
# Show simulation time
simulationTimeModelAnonGauss <- Ie-Is
print(Ie-Is)

#### Model B
# Model specific parameters
mu_ModelB    = function( x ){ sin( 8 * pi * x ) * exp( -3 * x ) }
noise_ModelB = function( N,
                         x,
                         sigma = function(x){ rep(1, length(x)) }, df){
                          randNumber=function(n){ (rchisq(n, df=df)-df)/sqrt(2*df)}
                          GaussDensitySumNoise( N=N, x=x, sigma=sigma, randNumber=randNumber ) }
sigma_ModelB = function(x){ ((1-x-0.4)^2+1)/6 }
dfVec        = c(1,2,7,10,30,60)

# initialize the matrix for the covering rates
covRate_ModelBnonGauss_tGKF <- covRate_ModelBnonGauss_boot <- covRate_ModelBnonGauss_mboot <- covRate_ModelBnonGauss_mrboot <- matrix( NA, length(NVec), length(dfVec) )
colnames(covRate_ModelBnonGauss_tGKF) <- colnames(covRate_ModelBnonGauss_boot) <- colnames(covRate_ModelBnonGauss_mboot) <- colnames(covRate_ModelBnonGauss_mrboot)  <- dfVec
rownames(covRate_ModelBnonGauss_tGKF) <- rownames(covRate_ModelBnonGauss_boot) <- rownames(covRate_ModelBnonGauss_mboot)  <- rownames(covRate_ModelBnonGauss_mrboot)  <- NVec

# Simulation of the covering rates
Is <-Sys.time()
for(k in 1:length(dfVec)){
  covRate_ModelBnonGauss_tGKF[,k] = t(covRate_simulation( Msim, scenario="SpN", method=methodVec[1], param_method=param_method[[1]], level=level, N=NVec, mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=0, x=x, df=dfVec[k] ))
}
Ie <- Sys.time()
# Show simulation time
simulationTimeModelBnonGauss_tGKF <- Ie-Is
print(Ie-Is)

Is <-Sys.time()
for(k in 1:length(dfVec)){
  covRate_ModelBnonGauss_boot[,k] = t(covRate_simulation( Msim, scenario="SpN", method=methodVec[2], param_method=param_method[[2]], level=level, N=NVec, mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=0, x=x, df=dfVec[k] ))
}
Ie <- Sys.time()
# Show simulation time
simulationTimeModelBnonGauss_boot <- Ie-Is
print(Ie-Is)

Is <-Sys.time()
for(k in 1:length(dfVec)){
  covRate_ModelBnonGauss_mboot[,k] = t(covRate_simulation( Msim, scenario="SpN", method=methodVec[3], param_method=list(1=param_method[[3]]), level=level, N=NVec, mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=0, x=x, df=dfVec[k] ))
}
Ie <- Sys.time()
# Show simulation time
simulationTimeModelBnonGauss_mboot <- Ie-Is
print(Ie-Is)

Is <-Sys.time()
for(k in 1:length(dfVec)){
  covRate_ModelBnonGauss_mrboot[,k] = t(covRate_simulation( Msim, scenario="SpN", method=methodVec[4], param_method=param_method[4], level=level, N=NVec, mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=0, x=x, df=dfVec[k] ))
}
Ie <- Sys.time()
# Show simulation time
simulationTimeModelBnonGauss_mrboot <- Ie-Is
print(Ie-Is)

#### save the simulation
rm(k, i, Ie, Is)
save.image(paste(path_sim, date,"_SimulationSmoothNonGauss.RData"))
rm( covRate_ModelBnonGauss_tGKF, covRate_ModelBnonGauss_boot, covRate_ModelAnonGauss )

############ Plot the Results and error processes
# date = "11_27_2018"
# load(paste(path_sim, date,"_SimulationSmoothNonGauss.RData"))
pdfname <- paste( path_pics, "ResultsSimulationNonGaussModelA.pdf", sep = "" )
pdf( pdfname, title = pdfname, width = 1.1 * 10, height = 1.2 * 4 )
  par( mfrow = c( 1, 2 ), mar = c( 4.1, 5.1, 2.1, 2.1 ) )
  Nmethods = length( param_method )
  set.seed(5)
  # Plot sample paths of noise of Model A
  x <- seq( 0, 1, length.out=100 )
  y <- noise_ModelA(15, x=x, sigma=function(x){ rep(1, length(x)) })

  plot( NULL, xlim = c(0,1), ylim = range(y), ylab="", xlab=" ", cex.axis=caxis, cex.main=cmain, cex.lab=clab )
  matlines(x, y)


  ###### Plot Results of simulation
  # Plot Results of Model A
  colVec <- c("red", "darkmagenta", "cyan3","blue")
  cov.data = covRate_ModelAnonGauss
  model    = "A"
  ylow     = 90
  Nlength = length(NVec)
  axisPoints = 10*(1:Nlength)
  plot( NULL, xlim = range(axisPoints), ylim = c(ylow,100), ylab="Covering Rate [%]", xlab="Number of Samples [N]", xaxt='n', cex.lab=clab, cex.axis=caxis )
  axis(1, at=axisPoints, labels=NVec, cex.axis=caxis)

  ### nominal level
  abline(h=level*100, col="black", lty=1)
  ### 2 times standard error of nominal level
  abline( h = 100*(level+c(-2,2)*sqrt(level*(1-level)/Msim)), col="black", lty=2)

  ### Plot the simulation results
  for(k in 1:Nmethods){
    lines( axisPoints, 100*cov.data[,k], col=colVec[k], pch=k  )
    points( axisPoints, 100*cov.data[,k], col=colVec[k], pch=k  )
  }
  legend("bottomright", legend=c("Boots-t", "tGKF", "rMult-t", "gMult-t"), col=colVec[c(2,1,4,3)], pch=c(2,1,4,3), lty=rep(1,Nmethods), cex=1.5)
dev.off()

## Plot model B
pdfname <- paste( path_pics, "ResultsSimulationNonGaussModelB.pdf", seq = "" )
pdf( pdfname, title = pdfname, width = 1.1 * 10, height = 1.2 * 5 )
par( mfrow = c( 2, 3 ), mar = c( 2.6, 5.1, 3.1, 2.1 ) )
     Nmethods = length( param_method )
     axisPoints = 10*(1:Nlength)
     set.seed(5)
    # Plot sample paths of noise of Model B
    for( k in c(1,3,6) ){
      y <- noise_ModelB( N=15, x=x, sigma=function(x){ rep(1, length(x)) }, df=dfVec[k] )
      plot( NULL, xlim = c(0,1), ylim = c(-3.5,3.5), ylab="", xlab="", main=bquote(nu == .(dfVec[k])), cex.axis=caxis, cex.main=cmain, cex.lab=clab )
      matlines(x,y)
    }

    ###### Plot Results of simulation
    par(mar=c(4.1,5.1,2.1,2.1) )

    at.points = 10*1:Nlength
    for(l in 1:3){
          if( l == 1 ){
            cov.data = 100*covRate_ModelBnonGauss_tGKF
            colFun <- colorRamp(c("rosybrown1","darkred"), bias = 1, space = "rgb", interpolate = "linear")
            colorVec = colFun(seq(0,1,length.out=6))
            colorVec = apply(colorVec, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
            colorVec[4:5] = c("red", "red3")
            ylims = c(75,100)
            title = "tGKF"
            place.legend = "bottomright"
          }else if( l == 2 ){
            cov.data = 100*covRate_ModelBnonGauss_mrboot
            colFun <- colorRamp(c("lightblue1","darkblue"), bias = 1, space = "rgb", interpolate = "linear")
            colorVec = colFun(seq(0,1,length.out=6))
            colorVec = apply(colorVec, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
            colorVec[4:5] = c("royalblue1","blue")
            ylims = c(75,100)
            title = "rMult-t"
            place.legend = "bottomright"
          }else{
            cov.data = 100*covRate_ModelBnonGauss_boot
            colFun <- colorRamp(c("plum1","darkmagenta"), bias = 1, space = "rgb", interpolate = "linear")
            colorVec = colFun(seq(0,1,length.out=6))
            colorVec = apply(colorVec, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
            ylims = c(93.5,100)
            title = "Bootstrap-t"
            place.legend = "topright"
          }

          plot( NULL, xlim = range(axisPoints), ylim = ylims, ylab="Covering Rate [%]", xlab="Number of Samples [N]", xaxt='n', cex.lab=clab, cex.axis=caxis, main= title, cex.main=cmain )
          ### nominal level
          abline(h=95, col="black", lty=1)
          ### 2 times standard error of nominal level
          abline( h = 100*(level+c(-2,2)*sqrt((1-level)*level/Msim)), col="black", lty=2)

          for(nn in 1:length(dfVec) ){
            lines( at.points, cov.data[,nn], col=colorVec[nn], pch=nn  )
            points( at.points, cov.data[,nn], col=colorVec[nn], pch=nn  )
          }
          axis(1, at=at.points, labels=NVec, cex.axis=caxis)
          legend(place.legend, legend=c(expression(nu), dfVec), col=c("white", colorVec), pch=0:length(dfVec), lty=rep(1,3), cex=1.2)
          rm(cov.data, place.legend, ylims, title, colorVec)
    }
dev.off()



## Clear Workspace
rm(ylow, Nlength, width, height, model, methodVec, l, at.points, axisPoints, colVec)
rm( covRate_ModelAnonGauss, covRate_ModelBnonGauss_tGKF, covRate_ModelBnonGauss_boot, i, Ie, Is, covRate_ModelBnonGauss_boot, param_method, covRate_ModelBnonGauss_mrboot, covRate_ModelBnonGauss_mboot )
rm( x, level, Msim, NVec, sd_ObsNoise, methodVec, date, y )
rm( mu_ModelA, mu_ModelB, noise_ModelA,noise_ModelB, sigma_ModelA, sigma_ModelB, sd_ObsNoise, dfVec )
rm( simulationTimeModelBnonGauss_mrboot, simulationTimeModelBnonGauss_mboot , simulationTimeModelBnonGauss_boot,
    simulationTimeModelAnonGauss, simulationTimeModelBnonGauss, pdfname  )

################################################################################################
##                                                                                            ##
##   Simulations and tables for Average Width and Variance of different Quantile Estimators   ##
##                                                                                            ##
################################################################################################
date = "11_14_2018"
# Compute the T statistic
Tstat <- function(X, true.mu){
  sqrt(dim(X)[2])*(rowMeans(X)-true.mu)/sqrt(matrixStats::rowVars(X))
}
################## Simulation Parameter
MonteN        = 1e5
x             = seq(0,1, length.out=200)
NSampVec      = length(NVec)
Msim          = 1e3
################## Gaussian case
# Parametsers
noise    = GaussDensitySumNoise

############# Monte Carlo Simulation for the Correct Quantiles
Is <- Sys.time()
MonteTgauss <- matrix( NA, MonteN, NSampVec )
for(k in 1:MonteN){
  y = noise( N=max(NVec), x=x )
  for(l in 1:NSampVec){
    MonteTgauss[k,l] <- max(abs(Tstat( X=y[,1:NVec[l]], true.mu=rep(0, length(x)) )))
  }
}
Ie <- Sys.time()
Ie-Is
rm(l,k)

true.quantgauss <- apply(MonteTgauss, 2, function(vec) quantile(vec,level, type=8))
print(true.quantgauss)

hatquantilesgauss <- array( NA, dim=c( Msim, NSampVec, 9) )

Is<-Sys.time()
# Simulate the Distribution of the Quantiles
for( m in 1:Msim ){
  # Simulate the data
  y = noise( N=max(NVec), x=x )
  for( n in 1:NSampVec ){
    # Generate data samples
    data = y[,sample( 1:max(NVec), size = NVec[n] )]
    # Estimate the LKCs
    hatLKC = c(1, LKC_estim_direct( data ) )
    # Compute the Quantiles of SCBs
    hatquantilesgauss[m,n,1]  <- GKFquantileApprox( alpha = (1-level)/2, LKC=hatLKC, field="t", df=NVec[n] )
    hatquantilesgauss[m,n,2]  <- GKFquantileApprox( alpha = (1-level)/2, LKC=hatLKC, field="Gauss" )
    degras <- SCBmeanfd::scb.mean( x, t(data), bandwidth=diff(x)[1]/3, level=level, scbtype="normal", nboot=5e3 )
    hatquantilesgauss[m,n,3]  <- degras$qnorm
    hatquantilesgauss[m,n,4]  <- NonParametricBootstrap( A=data, params=list(Mboots=5e3, alpha=1-level, method="regular") )$q
    hatquantilesgauss[m,n,5]  <- NonParametricBootstrap( A=data, params=list(Mboots=5e3, alpha=1-level, method="t") )$q
    hatquantilesgauss[m,n,6]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="regular", weights="gauss") )$q
    hatquantilesgauss[m,n,7]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="t", weights="gauss") )$q
    hatquantilesgauss[m,n,8]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="regular", weights="rademacher") )$q
    hatquantilesgauss[m,n,9]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="t", weights="rademacher") )$q
    hatquantilesgauss[m,n,10]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="regular", weights="mammen") )$q
    hatquantilesgauss[m,n,11]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="t", weights="mammen") )$q
  }
}
Ie <- Sys.time()
SimulationTimeQuant = Ie-Is
rm(n,m, Ie, Is, y, data)

################## non Gaussian case
# Parameters
noise    = function(N, x, sigma=function(x){ rep(1, length(x)) }, df=7){
  randNumber=function(n){ (rchisq(n, df=df)-df)/sqrt(2*df)}
  GaussDensitySumNoise( N=N, x=x, sigma=sigma, randNumber=randNumber ) }

############# Monte Carlo Simulation for the Correct Quantiles
Is <- Sys.time()
MonteTnongauss <- matrix( NA, MonteN, NSampVec )
for(k in 1:MonteN){
  y = noise( N=max(NVec), x=x )
  for(l in 1:NSampVec){
    MonteTnongauss[k,l] <- max(abs(Tstat( X=y[,1:NVec[l]], true.mu=rep(0, length(x)) )))
  }
}
Ie <- Sys.time()
Ie-Is
rm(l,k)

true.quantnongauss <- apply(MonteTnongauss, 2, function(vec) quantile(vec,level, type=8))
print(true.quantnongauss)

hatquantilesnongauss <- array( NA, dim=c( Msim, NSampVec, 9) )

Is<-Sys.time()
# Simulate the Distribution of the Quantiles
for( m in 1:Msim ){
  # Simulate the data
  y = noise( N=max(NVec), x=x )
  for( n in 1:NSampVec ){
    # Generate data samples
    data = y[,sample( 1:max(NVec), size = NVec[n] )]
    # Estimate the LKCs
    # hatLKC = c(1, LKC_estim_direct( data ) )
    # Compute the Quantiles of SCBs
    # hatquantilesnongauss[m,n,1]  <- GKFquantileApprox( alpha = (1-level)/2, LKC=hatLKC, field="t", df=NVec[n] )
    # hatquantilesnongauss[m,n,2]  <- GKFquantileApprox( alpha = (1-level)/2, LKC=hatLKC, field="Gauss" )
    # degras <- SCBmeanfd::scb.mean( x, t(data), bandwidth=diff(x)[1]/3, level=level, scbtype="normal", nboot=5e3 )
    # hatquantilesnongauss[m,n,3]  <- degras$qnorm
    # hatquantilesnongauss[m,n,4]  <- NonParametricBootstrap( A=data, params=list(Mboots=5e3, alpha=1-level, method="regular") )$q
    # hatquantilesnongauss[m,n,5]  <- NonParametricBootstrap( A=data, params=list(Mboots=5e3, alpha=1-level, method="t") )$q
    # hatquantilesnongauss[m,n,6]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="regular") )$q
    # hatquantilesnongauss[m,n,7]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="t") )$q
    # hatquantilesnongauss[m,n,8]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="regular", weights="rademacher") )$q
    # hatquantilesnongauss[m,n,9]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="t", weights="rademacher") )$q
    hatquantilesnongauss[m,n,10]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="regular", weights="mammen") )$q
    hatquantilesnongauss[m,n,11]  <- MultiplierBootstrap(R=data, params=list(Mboots=5e3, alpha=1-level, method="t", weights="mammen") )$q
  }
}
Ie <- Sys.time()
SimulationTimeQuantnon = Ie-Is
rm(n,m, Ie, Is, y, data)

## Save Workspace
save.image(paste(path_sim, date,"_SimulationQuantiles.RData",sep="") )
rm(hatquantilesnongauss, hatquantilesgauss, MonteTnongauss, true.quantnongauss, MonteTgauss, true.quantgauss )

#################################### Produce Latex table from Simulation ####################################
date = "08_17_2018"
load(paste(path_sim, date,"_SimulationQuantiles.RData",sep="") )

meansGauss <- apply(hatquantilesgauss, 2:3, mean)
varsGauss <- apply(hatquantilesgauss, 2:3, function(y) sqrt(var(y)) )

tableGauss = cbind( as.character(round(true.quantgauss,3) ), matrix(paste("$", as.character(round(meansGauss,3)),"\\pm", as.character(round(2*varsGauss/sqrt(NVec),3)),"$"),dim(meansGauss)[1],dim(meansGauss)[2]) )

colnames(tableGauss) <- c("\\textbf{true}", "\\textbf{tGKF}", "\\textbf{GKF}", "\\textbf{Degras}", "\\textbf{Boots}", "\\textbf{Boots-t}", "\\textbf{gMult}", "\\textbf{gMult-t}", "\\textbf{rMult}", "\\textbf{rMult-t}")
rownames(tableGauss) <- paste("$ \\mathbf{", as.character(NVec),"}$")
print(xtable( t(tableGauss) ), align="|l|cccccc|", sanitize.text.function=function(x){x} )

meansnonGauss <- apply(hatquantilesnongauss, 2:3, mean)
varsnonGauss <- apply(hatquantilesnongauss, 2:3, function(y) sqrt(var(y)) )
tablenonGauss = cbind( as.character(round(true.quantnongauss,3) ), matrix(paste("$", as.character(round(meansnonGauss,3)),"\\pm", as.character(round(2*varsnonGauss/sqrt(NVec),3)),"$"),dim(meansnonGauss)[1],dim(meansnonGauss)[2]) )
colnames(tablenonGauss) <- c("\\textbf{true}", "\\textbf{tGKF}", "\\textbf{GKF}", "\\textbf{Degras}", "\\textbf{Boots}", "\\textbf{Boots-t}", "\\textbf{gMult}", "\\textbf{gMult-t}", "\\textbf{rMult}", "\\textbf{rMult-t}")
rownames(tablenonGauss) <- paste("$ \\mathbf{", as.character(NVec),"}$")
print(xtable( t(tablenonGauss) ), align="|l|cccccc|", sanitize.text.function=function(x){x} )


################################################################################################
##                                                                                            ##
##               Simulations and Plots for Influence of Observation Noise                     ##
##                                                                                            ##
################################################################################################
############################### Observation and bandwidth dependence ###########################
############ General parameters
x           = seq( 0, 1, length.out=100 )
level       = 0.95
Msim        = 2.5e3
NVec        = c( 10, 20, 30, 50, 100, 150 )
N           = length(NVec)
date        = "12_05_2018"
obsSDvec    = c(0.05, 0.1, 0.2)
NobsSd      = length(obsSDvec)
BWvec       = c(0.02, 0.03, 0.05, 0.1)
xeval       = seq( 0, 1, length.out=200 )
############ Model A
# Model specific parameters
mu_ModelA    = function(x){ sin(8*pi*x) * exp(-3*x) }
noise_ModelA = BernsteinSumNoise
sigma_ModelA = function(x){ ((1-x-0.4)^2+1)/6 }

# initialize the matrix for the covering rates
covRate_ModelANonSmooth <- array( NA, dim = c(length(obsSDvec), length(NVec), length(BWvec)) )

Is <- Sys.time()
## Simulation of the covering rates
for( k in 1:length(BWvec) ){
  # precompute the smoothing matrix
  param_scenario <- list()
  param_scenario$xeval         = xeval
  param_scenario$SmoothWeights = locpol::locLinWeightsC( x=x, xeval=param_scenario$xeval , bw=BWvec[k], kernel=locpol::gaussK)$locWeig
  for( j in 1:length(NVec) ){for(i in 1:length(obsSDvec)){
    param_scenario$SmoothWeights
    tmp = covRate_simulation( Msim, scenario="SpN", param_scenario = param_scenario, method="tGKF", level=level, N=NVec[j], mu=mu_ModelA, noise=noise_ModelA, sigma=sigma_ModelA, sd_ObsNoise=obsSDvec[i], x=x )
    covRate_ModelANonSmooth[i,j,k] = tmp$covRateSmoothed
  }}
}
Ie <- Sys.time()
simulationTimeModelAnonSmooth = Ie-Is


############ Model B
# Model specific parameters
mu_ModelB    = function(x){ sin(8*pi*x) * exp(-3*x) }
noise_ModelB = GaussDensitySumNoise
sigma_ModelB = function(x){ ((1-x-0.4)^2+1)/6 }

# initialize the matrix for the covering rates
covRate_ModelBNonSmooth <- array( NA, dim = c(length(obsSDvec), length(NVec), length(BWvec)) )

Is <- Sys.time()
# Simulation of the covering rates
for( k in 1:length(BWvec) ){
  # precompute the smoothing matrix
  param_scenario <- list()
  param_scenario$xeval         = xeval
  param_scenario$SmoothWeights = locpol::locLinWeightsC( x=x, xeval=param_scenario$xeval, bw=BWvec[k], kernel=locpol::gaussK)$locWeig
  for( j in 1:length(NVec) ){for(i in 1:length(obsSDvec)){
    param_scenario$SmoothWeights
    covRate_ModelBNonSmooth[i,j,k] = covRate_simulation( Msim, scenario="SpN", param_scenario = param_scenario, method="tGKF", level=level, N=NVec[i], mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=obsSDvec[i], x=x )$covRateSmoothed
  }}
}
Ie <- Sys.time()
simulationTimeModelBnonSmooth = Ie-Is

#### save the simulation
rm(k, i, j, tmp, Ie, Is)
save.image( paste(path_sim, date,"_SimulationNonSmoothGauss.RData") )


############################### Scale Space Simulations ###########################
# initialize the matrix for the covering rates
covRate_ModelBScale <- array( NA, dim = c(length(NVec), length(obsSDvec)) )

# precompute the smoothing matrix
param_scenario        = list()
param_scenario$bw     = exp( seq( log(.02), log(0.1),length.out=100 ) )
kernel                = locpol::gaussK
param_scenario$kernel = kernel
param_scenario$xeval  = xeval

# Simulation of the covering rates
Is <- Sys.time()
for( j in 1:length(NVec) ){for(i in 1:length(obsSDvec)){
  param_scenario$SmoothWeights
  covRate_ModelBScale[j,i] = covRate_simulation( Msim, scenario="Scale1d", param_scenario = param_scenario, method="tGKF", level=level, N=NVec[j], mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=obsSDvec[i], x=x )
}}
Ie <- Sys.time()
simulationTimeModelBScale = Ie-Is

#### save the simulation
rm(i, j)
save.image( paste(path_sim, date,"_SimulationNonSmoothGauss.RData") )


############ Plot Results
load(paste(path_sim, date,"_SimulationNonSmoothGauss.RData"))
pdfname <- "Pics/ResultsSimulationInfluenceObsNoiseModelA.pdf"
pdf(pdfname,  title=pdfname, width=1.1*10, height=1.1*5)
par(mfrow=c(2,NobsSd), mar=c(2.1,5.1,5.1,2.1) )
sub = 1:length(BWvec)
h_Vec=BWvec[sub]
#level=0.98
for(k in 1:NobsSd){
  if( k!=2 ){
    title = ''
  }else{
    title = "Model A plus Noise"
  }
  x <- seq( 0, 1, length.out=100 )
  y = FunctionalDataSample( N=15, x=x, mu=mu_ModelA, noise=noise_ModelA, sigma=sigma_ModelA, sd_ObsNoise=obsSDvec[k] )
  plot( NULL, xlim = c(0,1), ylim = c(-1.5,1.5), ylab="", xlab="", main=title, cex.axis=caxis, cex.main=cmain, cex.lab=clab )
  text(x=0.5, y=1.2, labels= bquote(sigma == .(obsSDvec[k])), cex=1.5 )
  matlines(x,y)
}
par(mar=c(5.1,5.1,2.1,2.1) )
for(k in 1:3){
  # if( k!=2 ){
  #   title = ''
  # }else{
  #   title = "Model 1 with observation Noise"
  # }
  axisPoints = 10*(1:N)
  plot( NULL, xlim = range(axisPoints), ylim = c(94,100), ylab="Covering Rate [%]", xlab="Number of Samples [N]", main=title, xaxt='n', cex.axis=caxis, cex.main=cmain, cex.lab=clab )
  axis(1, at=axisPoints, labels=NVec, cex.axis=caxis)

  ### nominal level
  abline(h=level*100, col="black", lty=1)
  ### 2 times standard error of nominal level
  abline( h = 100*(level+c(-2,2)*sqrt(level*(1-level)/Msim)), col="black", lty=2)

  ### Plot the data
  cov.data <- covRate_ModelANonSmooth[k,,sub]
  for( j in 1:dim(cov.data)[2] ){
    #   text(x=20, y=98.8, labels= bquote(sigma == .(Noise_sdVec[k])), cex=1.8 )
    lines(  axisPoints, cov.data[,j], col=j, pch=k, lty=1  )
    points( axisPoints, cov.data[,j], col=j, pch=j  )
  }
  if(k==1){
    legend("topright", legend=c("Bandwidth", as.character(h_Vec)), col=c("white",1:length(h_Vec)), pch=c(0:length(h_Vec)), lty=rep(1,3), cex=1.2)
  }
}
dev.off()

pdfname <- "Pics/ResultsSimulationInfluenceObsNoiseModelB.pdf"
pdf(pdfname,  title=pdfname, width=1.1*10, height=1.1*5)
par(mfrow=c(2,NobsSd), mar=c(2.1,5.1,5.1,2.1) )
sub = 1:length(BWvec)
h_Vec=BWvec[sub]
#level=0.98
for(k in 1:NobsSd){
  if( k!=2 ){
    title = ''
  }else{
    title = "Model B plus Noise"
  }
  x <- seq( 0, 1, length.out=100 )
  y = FunctionalDataSample( N=15, x=x, mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=obsSDvec[k] )
  plot( NULL, xlim = c(0,1), ylim = c(-1.5,1.5), ylab="", xlab="", main=title, cex.axis=caxis, cex.main=cmain, cex.lab=clab )
  text(x=0.5, y=1.2, labels= bquote(sigma == .(obsSDvec[k])), cex=1.5 )
  matlines(x,y)
}
par(mar=c(5.1,5.1,2.1,2.1) )
for(k in 1:3){
  # if( k!=2 ){
  #   title = ''
  # }else{
  #   title = "Model 1 with observation Noise"
  # }
  axisPoints = 10*(1:N)
  plot( NULL, xlim = range(axisPoints), ylim = c(94,100), ylab="Covering Rate [%]", xlab="Number of Samples [N]", main=title, xaxt='n', cex.axis=caxis, cex.main=cmain, cex.lab=clab )
  axis(1, at=axisPoints, labels=NVec, cex.axis=caxis)

  ### nominal level
  abline(h=level*100, col="black", lty=1)
  ### 2 times standard error of nominal level
  abline( h = 100*(level+c(-2,2)*sqrt(level*(1-level)/Msim)), col="black", lty=2)

  ### Plot the data
  cov.data <- covRate_ModelBNonSmooth[k,,sub]
  for( j in 1:dim(cov.data)[2] ){
    #   text(x=20, y=98.8, labels= bquote(sigma == .(Noise_sdVec[k])), cex=1.8 )
    lines(  axisPoints, cov.data[,j], col=j, pch=k, lty=1  )
    points( axisPoints, cov.data[,j], col=j, pch=j  )
  }
  if(k==1){
    legend("topright", legend=c("Bandwidth", as.character(h_Vec)), col=c("white",1:length(h_Vec)), pch=c(0:length(h_Vec)), lty=rep(1,3), cex=1.2)
  }
}
dev.off()

# Plot Scale Space Results
pdfname <- "Pics/ResultsSimulationScaleGaussianB.pdf"
pdf(pdfname,  title=pdfname, width=1.1*10, height=1*3)
par(mfrow=c(1,3))
par(mar=c(2.1,5.1,3.1,2.1) )
# Plot Scale space noise amples Model B
set.seed(pi)
x = seq(0,1, length.out=100)

## Simulate noise and plot
Y   = FunctionalDataSample( N=2, x=x, mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=obsSDvec[1] )
Y   = scaleField( Y, x=x, xeval=param_scenario$xeval, h=param_scenario$bw )
zlim=range( c(as.vector(Y[,,2]), as.vector(Y[,,1]) ))
fields::image.plot(Y[,,2], cex.axis=caxis, cex.main=cmain, zlim=zlim, main = bquote(sigma == .(obsSDvec[1])))

Y   = FunctionalDataSample( N=2, x=x, mu=mu_ModelB, noise=noise_ModelB, sigma=sigma_ModelB, sd_ObsNoise=obsSDvec[3] )
Y   = scaleField( Y, x=x, xeval=param_scenario$xeval, h=param_scenario$bw )
zlim=range( c(as.vector(Y[,,2]), as.vector(Y[,,1]) ))
fields::image.plot(Y[,,2], cex.axis=caxis, cex.main=cmain, zlim=zlim, main = bquote(sigma == .(obsSDvec[3])))

# Plot Results
colVec <- c("darkorchid", "blue", "cyan3")
cov.data = covRate_ModelBScale
model    = "B"

axisPoints = 10*(1:N)
plot( NULL, xlim = range(axisPoints), ylim = c(94,96), ylab="Covering Rate [%]", xlab="Number of Samples [N]", xaxt='n', cex.lab=clab, cex.axis=caxis )
axis(1, at=axisPoints, labels=NVec, cex.axis=caxis)

### nominal level
abline(h=level*100, col="black", lty=1)
### 2 times standard error of nominal level
abline( h = 100*(level+c(-2,2)*sqrt(level*(1-level)/Msim)), col="black", lty=2)

### Plot the data
for(k in 1:3){
  lines(  axisPoints, 100*cov.data[,k], col=colVec[k], pch=k  )
  points( axisPoints, 100*cov.data[,k], col=colVec[k], pch=k  )
}

legend("bottomleft", legend=c(expression(sigma), obsSDvec), col=c("white", colVec), pch=0:length(obsSDvec), lty=rep(1,3), cex=1.2)

dev.off()





rm(covRate_ModelANonSmooth, covRate_ModelBNonSmooth, covRate_ModelBScale)


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
# Msim      = 5e3
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
#     Msim, N1 = N[1,i], N2 = N[2,i], level = level, x = x,
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
#     Msim, N1 = N[1,i], N2 = N[2,i], level = level, x = x,
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
#     Msim, N1 = N[1,i], N2 = N[2,i], level = level, x = x,
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
#     Msim, N1 = N[1,i], N2 = N[2,i], level = level, x = x,
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
#   abline( h = 100*0.95+100*c(1,-1)*2*sqrt(0.95*0.05/Msim), col="black", lty=2)
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
#   abline( h = 95+100*c(1,-1)*2*sqrt(0.95*0.05/Msim), col="black", lty=2)
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
#   abline( h = 95+100*c(1,-1)*2*sqrt(0.95*0.05/Msim), col="black", lty=2)
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
