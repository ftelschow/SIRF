# Clear workspace
rm( list = ls() )

# Required packages
require(RFT)
require(SampleFields)
require(tidyverse)
require(gaussplotR)
require(fields)
require(reshape2)

#
path_wd   <- "~/Seafile/Code/Rpackages/SIRF/Code_Articles/FunctionalDeltaResiduals/"
path_pics <- "~/Seafile/Projects/2019_DeltaResiduals/Article/Figures/"
# Set the working path
setwd(path_wd)

# Get the simulation file from source
#source(paste(path_wd, "Workspaces/Variance_Simulation.Rdata", sep = ""))

date  = "2021_11_15"
model = "ModelC" #  "ModelA" #  "ModelC" #
transformation.l ="cohensd" #  "skewness" #

if( transformation.l == "skewness" ){
  ylow = 0.75
  yup  = 1
}else{
  ylow = 0.85
  yup  = 1
}

#-------------------------------------------------------------------------------
sLegend <- 25
sText   <- 25
Sylab = 1.5
Sxlab = 1.5
sLine = 1.5
sPch = 4

theme1 <- theme(legend.position = "none",
                plot.title   = element_text(face = "bold", size = sText),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                axis.title.y = element_text(color = "black", size = sText, face = "plain"))

theme2 <- theme(legend.title = element_blank(),
                legend.position = c(0.82, 0.3),
                legend.key.size = unit(1.1, 'cm'),
                legend.text = element_text(size = sLegend),
                plot.title   = element_text(face = "bold", size = sText),
                axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                axis.title.y = element_text(color = "black", size = sText, face = "plain"))



#-------------------------------------------------------------------------------
# FALSE / FALSE plot
bias.l   = FALSE
se.est.l = FALSE

title = "No Bias correction, Gaussian variance"

# Collecting the data
load( file = paste(path_wd, "Workspaces/",
                date, "_",
                model, "_",
                transformation.l,
                "_bias_", bias.l,
                "_se.est_", se.est.l,
                "_full.rdata",
                sep = "" ) )

# Plot the simulation results
pngname <- paste( path_pics,
                  model, "_",
                  transformation.l,
                  "_bias_", bias.l,
                  "_se.est_", se.est.l,
                  ".png",
                  sep = "" )
png( pngname, width = 550, height = 450 )
target <- sqrt( lvl * ( 1 - lvl ) / Msim ) * c(-qnorm(lvl/2 + 1/2), 0, qnorm(lvl/2 + 1/2)) + lvl
xLab <- rownames(covRates)
yLab <- colnames(covRates)
covs <- as_tibble(covRates, rownames = "N") %>%
  melt(., id.vars = "N", variable = "Method", value.name = "CovRate") %>%
  as_tibble() %>% mutate_if(is.character, as.numeric)

# Plot the Covering Rates by Method
print( ggplot(covs, aes(N, CovRate, group = Method, col = Method)) +
         xlab( "Sample Size [N]" ) +
         ylab( "Covering Rate" ) +
         ggtitle( title ) +
         geom_hline( yintercept = target[2], size = sLine ) +
         geom_hline( yintercept = target[1], size = sLine, linetype = "twodash" ) +
         geom_hline( yintercept = target[3], size = sLine, linetype = "twodash" ) +
         geom_point(size = sPch) + geom_line(size = sLine) +
         theme2 +
         coord_cartesian(ylim = c(ylow, yup ))
)
dev.off()

#-------------------------------------------------------------------------------
# FALSE / TRUE plot
bias.l   = FALSE
se.est.l = TRUE

title = "No Bias correction, Estimated variance"

# Collecting the data
load( file = paste(path_wd, "Workspaces/",
                   date, "_",
                   model, "_",
                   transformation.l,
                   "_bias_", bias.l,
                   "_se.est_", se.est.l,
                   "_full.rdata",
                   sep = "" ) )

# Plot the simulation results
pngname <- paste( path_pics,
                  model, "_",
                  transformation.l,
                  "_bias_", bias.l,
                  "_se.est_", se.est.l,
                  ".png",
                  sep = "" )
png( pngname, width = 550, height = 450 )
target <- sqrt( lvl * ( 1 - lvl ) / Msim ) * c(-qnorm(lvl/2 + 1/2), 0, qnorm(lvl/2 + 1/2)) + lvl
xLab <- rownames(covRates)
yLab <- colnames(covRates)
covs <- as_tibble(covRates, rownames = "N") %>%
  melt(., id.vars = "N", variable = "Method", value.name = "CovRate") %>%
  as_tibble() %>% mutate_if(is.character, as.numeric)

# Plot the Covering Rates by Method
if( model == "ModelC" & transformation.l == "cohensd" ){
  print( ggplot(covs, aes(N, CovRate, group = Method, col = Method)) +
           xlab( "Sample Size [N]" ) +
           ylab( " " ) +
           ggtitle( title ) +
           geom_hline( yintercept = target[2], size = sLine ) +
           geom_hline( yintercept = target[1], size = sLine, linetype = "twodash" ) +
           geom_hline( yintercept = target[3], size = sLine, linetype = "twodash" ) +
           geom_point(size = sPch) + geom_line(size = sLine) +
           theme2 +
           coord_cartesian(ylim = c(ylow, yup ))
  )
}else{
  print( ggplot(covs, aes(N, CovRate, group = Method, col = Method)) +
           xlab( "Sample Size [N]" ) +
           ylab( " " ) +
           ggtitle( title ) +
           geom_hline( yintercept = target[2], size = sLine ) +
           geom_hline( yintercept = target[1], size = sLine, linetype = "twodash" ) +
           geom_hline( yintercept = target[3], size = sLine, linetype = "twodash" ) +
           geom_point(size = sPch) + geom_line(size = sLine) +
           theme1 +
           coord_cartesian(ylim = c(ylow, yup ))
  )
}
dev.off()

#-------------------------------------------------------------------------------
# TRUE / FALSE plot
bias.l   = TRUE
se.est.l = FALSE

title = "Bias correction, Gaussian variance"

# Collecting the data
load( file = paste(path_wd, "Workspaces/",
                   date, "_",
                   model, "_",
                   transformation.l,
                   "_bias_", bias.l,
                   "_se.est_", se.est.l,
                   "_full.rdata",
                   sep = "" ) )

# Plot the simulation results
pngname <- paste( path_pics,
                  model, "_",
                  transformation.l,
                  "_bias_", bias.l,
                  "_se.est_", se.est.l,
                  ".png",
                  sep = "" )
png( pngname, width = 550, height = 450 )
target <- sqrt( lvl * ( 1 - lvl ) / Msim ) * c(-qnorm(lvl/2 + 1/2), 0, qnorm(lvl/2 + 1/2)) + lvl
xLab <- rownames(covRates)
yLab <- colnames(covRates)
covs <- as_tibble(covRates, rownames = "N") %>%
  melt(., id.vars = "N", variable = "Method", value.name = "CovRate") %>%
  as_tibble() %>% mutate_if(is.character, as.numeric)

# Plot the Covering Rates by Method
print( ggplot(covs, aes(N, CovRate, group = Method, col = Method)) +
         xlab( "Sample Size [N]" ) +
         ylab( " " ) +
         ggtitle( title ) +
         geom_hline( yintercept = target[2], size = sLine ) +
         geom_hline( yintercept = target[1], size = sLine, linetype = "twodash" ) +
         geom_hline( yintercept = target[3], size = sLine, linetype = "twodash" ) +
         geom_point(size = sPch) + geom_line(size = sLine) +
         theme1 +
         coord_cartesian(ylim = c(ylow, yup ))
)
dev.off()

#-------------------------------------------------------------------------------
# TRUE / TRUE plot
bias.l   = TRUE
se.est.l = TRUE

title = "Bias correction, Estimated variance"

# Collecting the data
load( file = paste(path_wd, "Workspaces/",
                   date, "_",
                   model, "_",
                   transformation.l,
                   "_bias_", bias.l,
                   "_se.est_", se.est.l,
                   "_full.rdata",
                   sep = "" ) )

# Plot the simulation results
pngname <- paste( path_pics,
                  model, "_",
                  transformation.l,
                  "_bias_", bias.l,
                  "_se.est_", se.est.l,
                  ".png",
                  sep = "" )
png( pngname, width = 550, height = 450 )
target <- sqrt( lvl * ( 1 - lvl ) / Msim ) * c(-qnorm(lvl/2 + 1/2), 0, qnorm(lvl/2 + 1/2)) + lvl
xLab <- rownames(covRates)
yLab <- colnames(covRates)
covs <- as_tibble(covRates, rownames = "N") %>%
  melt(., id.vars = "N", variable = "Method", value.name = "CovRate") %>%
  as_tibble() %>% mutate_if(is.character, as.numeric)

# Plot the Covering Rates by Method
print( ggplot(covs, aes(N, CovRate, group = Method, col = Method)) +
         xlab( "Sample Size [N]" ) +
         ylab( " " ) +
         ggtitle( title ) +
         geom_hline( yintercept = target[2], size = sLine ) +
         geom_hline( yintercept = target[1], size = sLine, linetype = "twodash" ) +
         geom_hline( yintercept = target[3], size = sLine, linetype = "twodash" ) +
         geom_point(size = sPch) + geom_line(size = sLine) +
         theme1 +
         coord_cartesian(ylim = c(ylow, yup ))
)
dev.off()

