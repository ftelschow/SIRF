# Clear workspace
rm( list = ls() )

# Required packages
require(RFT)
require(SampleFields)
require(tidyverse)
require(gaussplotR)
require(fields)

# Load the data
path_wd   <- "/home/fabian/Seafile/Code/Rpackages/SCBfun/DeltaResiduals/"
load(paste(path_wd, "Workspaces/Variance_Simulation.Rdata", sep = ""))
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

#-------------------------------------------------------------------------------
# Theme of the plots
#-------------------------------------------------------------------------------
sLegend <- 25
sText   <- 25
Sylab = 1.5 
Sxlab = 1.5
sLine = 2

theme1 <- theme(line = element_line(size = ),
                legend.position = "none",
                plot.title   = element_text(face = "bold", size = sText),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),  
                axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                axis.title.y = element_text(color = "black", size = sText, face = "plain"))

theme2 <- theme(line = element_line(size = ),
                legend.title = element_text(face = "bold", size = sText),
                legend.position = c(0.82, 0.3),
                legend.key.size = unit(1.1, 'cm'),
                legend.text = element_text(size = sLegend),
                plot.title   = element_text(face = "bold", size = sText),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),  
                axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                axis.title.y = element_text(color = "black", size = sText, face = "plain"))


# Cohen's d variance estimation plot
pngname <- paste( path_pics,
                  "CohensdVarianceEstimate.png",
                  sep = "" )
png( pngname, width = 550, height = 450 )
print( ggplot(relBias_d, aes(x, rel.bias, group = N, col = N)) +
         geom_line(size = sLine) +
         xlab( "Location" ) +
         ylab( "Relative Bias" ) +
         ggtitle( "Cohen's d: Variance Estimate" ) +
         theme2) + 
         coord_cartesian(ylim = c(-0.05, 0)) +
         labs(colour = "Sample Size") 
dev.off()

# Skewness variance estimation plot
pngname <- paste( path_pics,
                  "SkewnessVarianceEstimate.png",
                  sep = "" )
png( pngname, width = 550, height = 450 )
print( ggplot(relBias_s, aes(x, rel.bias, group = N, col = N)) +
         geom_line(size = sLine) +
         xlab( "Location" ) +
         ylab( "Relative Bias" ) +
         ggtitle( "Skewness: Variance Estimate" ) +
         theme1 + 
         coord_cartesian(ylim = c(-0.30, 0)) )
dev.off()
