# Clear workspace
rm( list = ls() )

# Required packages
require(SampleFields)
require(tidyverse)

# Paths
path_wd   <- "/home/fabian/Seafile/Code/Rpackages/SCBfun/DeltaResiduals/"
path_pics <- "/home/fabian/Seafile/Projects/2019_DeltaResiduals/Article/Figures/"

# Set the working path
setwd(path_wd)

# coordinates
x = seq(0, 1, length.out = 100 )
N = 10

# Theme for the plots
sLegend <- 25
sText   <- 25
sTitle  <- 35
sLineMean   <- 3
sLine   <- 1.5
theme1 <- theme(legend.position = "none",
                plot.title   = element_text(face = "bold", size = sTitle),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                axis.title.x = element_text(color = "black", size = sText , face = "plain"),
                axis.title.y = element_text(color = "black", size = sText , face = "plain"),
)

# Model A
mu_model    = function(x){sin(4 * pi * x) * exp(-3 * x)}
sigma_model = function(x){((1 - x - 0.4)^2 + 1) / 6}
noise_model = RandomNormalSum

Y1 <- SignalPlusNoise(N  = N,
                      x  = x,
                      mu = mu_model,
                      noise = noise_model,
                      sigma = sigma_model)

tmp = plot_RF(Y1, size.line = sLine)
tbl = tibble(x = x, y = mu_model(x))

# Plot the simulation results
pngname <- paste( path_pics,
                  "ModelA_Samples.png",
                  sep = "" )
png( pngname, width = 550, height = 450 )
tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) + theme1 +
  ggtitle( "Model A" )
dev.off()


# ModelB
mu_model    = function(x){ ( x-0.3 )^2 }
sigma_model = function(x){ (sin(3*pi*x) + 1.5) / 6 }

covf <- function(x, y) SampleFields::covf.nonst.matern(x,
                                                       y,
                                                       params = c( 1, 1 / 4, 0.4 ))
noise_model = function(N, x){ ArbCovProcess( N, x, covf = covf ) }

Y1 <- SignalPlusNoise(N  = N,
                      x  = x,
                      mu = mu_model,
                      noise = noise_model,
                      sigma = sigma_model)

tmp = plot_RF(Y1, size.line = sLine)
tbl = tibble(x = x, y = mu_model(x))

# Plot the simulation results
pngname <- paste( path_pics,
                  "ModelB_Samples.png",
                  sep = "" )
png( pngname, width = 550, height = 450 )
tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) + theme1 +
  ggtitle( "Model B" )
dev.off()


# ModelC
mu_model    = function(x){sin(4 * pi * x) * exp(-3 * x)}
sigma_model = function(x){(1.5 - x)/2}
noise_model = DegrasNonGaussProcess

Y1 <- SignalPlusNoise(N  = 20,
                      x  = x,
                      mu = mu_model,
                      noise = noise_model,
                      sigma = sigma_model)
tmp = plot_RF(Y1, size.line = sLine)
tbl = tibble( x = x, y = mu_model(x) )

# Plot the simulation results
pngname <- paste( path_pics,
                  "ModelC_Samples.png",
                  sep = "" )
png( pngname, width = 550, height = 450 )
tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) + theme1 +
  ggtitle( "Model C" )
dev.off()


# Model A with observation noise
mu_model    = function(x){sin(4 * pi * x) * exp(-3 * x)}
sigma_model = function(x){((1 - x - 0.4)^2 + 1) / 6}
noise_model = RandomNormalSum

x = seq(0, 1, length.out = 50)

Y1 <- SignalPlusNoise(N  = N,
                      x  = x,
                      mu = mu_model,
                      noise = noise_model,
                      sigma = sigma_model,
                      obs.noise = 0.05)

tmp = plot_RF(Y1, size.line = sLine)
tbl = tibble(x = x, y = mu_model(x))

# Plot the simulation results
pngname <- paste( path_pics,
                  "ModelA_Samples05.png",
                  sep = "" )
png( pngname, width = 550, height = 450 )
tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) + theme1
dev.off()

Y1 <- SignalPlusNoise(N  = N,
                      x  = x,
                      mu = mu_model,
                      noise = noise_model,
                      sigma = sigma_model,
                      obs.noise = 0.1)

tmp = plot_RF(Y1, size.line = sLine)
tbl = tibble(x = x, y = mu_model(x))

# Plot the simulation results
pngname <- paste( path_pics,
                  "ModelA_Samples1.png",
                  sep = "" )
png( pngname, width = 550, height = 450 )
tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) + theme1
dev.off()
