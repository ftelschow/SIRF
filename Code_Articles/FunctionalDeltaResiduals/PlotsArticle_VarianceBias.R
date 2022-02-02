# Clear workspace
rm( list = ls() )

# Required packages
require(RFT)
require(SampleFields)
require(tidyverse)
require(gaussplotR)
require(fields)

hvalue = 6
wvalue = 1.2*8

# Load the data
path_wd   <- "/home/fabian/Seafile/Code/Rpackages/SCBfun/DeltaResiduals/"
load(paste(path_wd, "Workspaces/Variance_Simulation.Rdata", sep = ""))

path_wd   <- "/home/fabian/Seafile/Code/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/"
load(paste(path_wd, "Workspaces/Variance_Simulation2.Rdata", sep = ""))
path_pics <- "/home/fabian/Seafile/Projects/2019_DeltaResiduals/Article/Figures/"

# Compute the relative biases
names <- as.character(x)

# make relative bias into a tibble object
relBias_d = (apply( sim_var_d, 1:2, mean ) - trueVar_d) / trueVar_d
rownames(relBias_d) <- Nvec
relBias_d <- as_tibble(relBias_d, rownames = "N") %>%
  rename_at(vars(starts_with('V')), ~ names) %>%
  pivot_longer(!N, names_to = "x", values_to = "rel.bias") %>%
  mutate( N = as_factor(N)) %>%
  mutate( x = as.numeric(x))

# make relative bias into a tibble object
relBias_s = (apply( sim_var_skew, 1:2, mean ) - trueVar_skew ) / trueVar_skew
rownames(relBias_s) <- Nvec
relBias_s <- as_tibble(relBias_s, rownames = "N") %>%
  rename_at(vars(starts_with('V')), ~ names) %>%
  pivot_longer(!N, names_to = "x", values_to = "rel.bias") %>%
  mutate( N = as_factor(N)) %>%
  mutate( x = as.numeric(x))

# make relative bias into a tibble object
relBias_sn = (apply( sim_var_sn, 1:2, mean ) - trueVar_sn ) / trueVar_sn
rownames(relBias_sn) <- Nvec
relBias_sn <- as_tibble(relBias_sn, rownames = "N") %>%
  rename_at(vars(starts_with('V')), ~ names) %>%
  pivot_longer(!N, names_to = "x", values_to = "rel.bias") %>%
  mutate( N = as_factor(N)) %>%
  mutate( x = as.numeric(x))

#-------------------------------------------------------------------------------
# Theme of the plots
#-------------------------------------------------------------------------------
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

theme1 <- theme(line = element_line(size = ),
                legend.position = "none",
                plot.title   = element_text(face = "bold", size = sTitle),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                axis.title.y = element_text(color = "black", size = sText, face = "plain"))

theme2 <- theme(line = element_line(size = ),
                legend.title = element_text(face = "bold", size = sLegend),
                legend.position = c(0.82, 0.32),
                legend.key.size = unit(1.1, 'cm'),
                legend.text = element_text(size = sLegend),
                plot.title   = element_text(face = "bold", size = sTitle),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                axis.title.y = element_text(color = "black", size = sText, face = "plain"))


# Cohen's d variance estimation plot
pngname <- paste( path_pics,
                  "CohensdVarianceEstimate.png",
                  sep = "" )
X11(width = wvalue, height = hvalue)
print( ggplot(relBias_d, aes(x, rel.bias, group = N, col = N)) +
         geom_line(size = sLine) +
         xlab( "Location" ) +
         ylab( "Relative Bias" ) +
         ggtitle( "Cohen's d" ) +
         theme2) +
         coord_cartesian(ylim = c(-0.05, 0)) +
         labs(colour = "Sample Size")
savePlot(filename = pngname)
dev.off()

# Skewness variance estimation plot
pngname <- paste( path_pics,
                  "SkewnessVarianceEstimate.png",
                  sep = "" )
X11(width = wvalue, height = hvalue)
print( ggplot(relBias_s, aes(x, rel.bias, group = N, col = N)) +
         geom_line(size = sLine) +
         xlab( "Location" ) +
         ylab( "Relative Bias" ) +
         ggtitle( "Skewness" ) +
         theme1 +
         coord_cartesian(ylim = c(-0.30, 0)) )
savePlot(filename = pngname)
dev.off()

# Skewness (normality) variance estimation plot
pngname <- paste( path_pics,
                  "SkewnessNormalityVarianceEstimate.png",
                  sep = "" )
X11(width = wvalue, height = hvalue)
print( ggplot(relBias_sn, aes(x, rel.bias, group = N, col = N)) +
         geom_line(size = sLine) +
         xlab( "Location" ) +
         ylab( "Relative Bias" ) +
         ggtitle( "Transformed Skewness" ) +
         theme1 +
         coord_cartesian(ylim = c(-0.30, 0)) )
savePlot(filename = pngname)
dev.off()
