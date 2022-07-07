# Clear workspace
rm( list = ls() )

# Required packages
require(SampleFields)
require(tidyverse)
require(latex2exp)

# Paths
path_wd   <- "/home/fabian/Seafile/Code/Rpackages/SCBfun/DeltaResiduals/"
path_pics <- "/home/fabian/Seafile/Projects/2019_DeltaResiduals/Article/Figures/"

# Set the working path
setwd(path_wd)
# coordinates
N = 10

#Our transformation function
scaleFUN <- function(x) sprintf("%.2f", x)

sLegend <- 25
sText   <- 30
sTitle  <- 35
sLegend <- 35
Sylab = 1.5
Sxlab = 1.5
sLine = 1.5
sPch = 4

hvalue = 7
wvalue = 1.3*hvalue

sLineMean = 3

theme2 <- theme(legend.position = c(0.55, 0.88),
                legend.title = element_blank(),
                legend.text  = element_text(color = "black", size = sLegend, face = "plain"),
                legend.direction = "horizontal",
                legend.key = element_rect(colour = "transparent", fill = "transparent"),
                legend.background = element_blank(),
                plot.title   = element_text(face = "bold", size = sTitle),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                axis.title.y = element_text(color = "black", size = sText, face = "plain"))

theme1 <- theme(legend.position = "none",
                legend.title = element_blank(),
                legend.text  = element_text(color = "black", size = sText, face = "plain"),
                legend.direction = "horizontal",
                legend.key = element_rect(colour = "transparent", fill = "transparent"),
                legend.background = element_blank(),
                plot.title   = element_text(face = "bold", size = sTitle),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                axis.title.y = element_text(color = "black", size = sText, face = "plain"))

x = seq(0,1, length.out=200)

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
pdfname <- paste( path_pics,
                  "ModelA_Samples.pdf",
                  sep = "" )
X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) +
        scale_y_continuous(labels=scaleFUN) + theme1 +
        ggtitle( "Model A" ) +
        xlab( "Spatial Index [s]" ) + ylab( "" ))
savePlot(filename = pdfname)
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
pdfname <- paste( path_pics,
                  "ModelB_Samples.pdf",
                  sep = "" )
X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) +
        scale_y_continuous(labels=scaleFUN) + theme1 +
        ggtitle( "Model B" ) +
        xlab( "Spatial Index [s]" ) + ylab( "" ))
savePlot(filename = pdfname)
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
pdfname <- paste( path_pics,
                  "ModelC_Samples.pdf",
                  sep = "" )
X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) +
        scale_y_continuous(labels=scaleFUN) + theme1 +
        xlab( "Spatial Index [s]" ) + ylab( "" ) +
  ggtitle( "Model C" ))
savePlot(filename = pdfname)
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
pdfname <- paste( path_pics,
                  "ModelA_Samples05.pdf",
                  sep = "" )

X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) +
        scale_y_continuous(labels=scaleFUN) +
        xlab( "Spatial Index [s]" ) + ylab( "" ) + theme1 +
        ggtitle( TeX(r'(\textbf{Model A} ($\sigma = 0.05$))')))
savePlot(filename = pdfname)
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
pdfname <- paste( path_pics,
                  "ModelA_Samples1.pdf",
                  sep = "" )
X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) +
        scale_y_continuous(labels=scaleFUN) +
        xlab( "Spatial Index [s]" ) + ylab( "" ) + theme1 +
        ggtitle( TeX(r'(\textbf{Model A} ($\sigma = 0.1$))')))
savePlot(filename = pdfname)
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
pdfname <- paste( path_pics,
                  "ModelB_Samples05.pdf",
                  sep = "" )
X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) +
        scale_y_continuous(labels=scaleFUN) +
        xlab( "Spatial Index [s]" ) + ylab( "" ) + theme1 +
        ggtitle( TeX(r'(\textbf{Model B} ($\sigma = 0.05$))')))
savePlot(filename = pdfname)
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
pdfname <- paste( path_pics,
                  "ModelB_Samples1.pdf",
                  sep = "" )
X11(width = wvalue, height = hvalue)
print(tmp + geom_line( data = tbl, aes(x = x, y = y), size = sLineMean ) +
        scale_y_continuous(labels=scaleFUN) +
        xlab( "Spatial Index [s]" ) + ylab( "" ) + theme1 +
        ggtitle( TeX(r'(\textbf{Model B} ($\sigma = 0.1$) )')))
savePlot(filename = pdfname)
dev.off()
