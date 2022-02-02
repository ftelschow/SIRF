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
sLegend <- 25
sText   <- 30
sTitle  <- 35
sLegend <- 25
Sylab = 1.5
Sxlab = 1.5
sLine = 1.5
sPch = 4
sLineMean   <- 3

wvalue = 1.3*7
hvalue = 7

# coordinates
x = seq(0, 1, length.out = 100 )
N = 10

# Theme for the plots
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
X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) + theme1 +
        ggtitle( "Model A" ))
savePlot(filename = pngname)
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
X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) + theme1 +
        ggtitle( "Model B" ))
savePlot(filename = pngname)
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
X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) + theme1 +
  ggtitle( "Model C" ))
savePlot(filename = pngname)
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

X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) + theme1 +
        ggtitle( TeX(r'(\textbf{Model A} ($\sigma = 0.05$))')))
savePlot(filename = pngname)
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
X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) + theme1 +
        ggtitle( TeX(r'(\textbf{Model A} ($\sigma = 0.1$))')))
savePlot(filename = pngname)
dev.off()


# Model B with observation noise
mu_model    = function(x){ ( x - 0.3 )^2 }
sigma_model = function(x){ (sin(3 * pi * x) + 1.5) / 6 }

covf <- function(x, y) SampleFields::covf.nonst.matern(x,
                                                       y,
                                                       params = c( 1, 1 / 4, 0.4 ))
noise_model = function(N, x){ ArbCovProcess( N, x, covf = covf ) }

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
                  "ModelB_Samples05.png",
                  sep = "" )
X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) + theme1 +
        ggtitle( TeX(r'(\textbf{Model B} ($\sigma = 0.05$))')))
savePlot(filename = pngname)
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
                  "ModelB_Samples1.png",
                  sep = "" )
X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) + theme1 +
        ggtitle( TeX(r'(\textbf{Model B} ($\sigma = 0.1$) )')))
savePlot(filename = pngname)
dev.off()
