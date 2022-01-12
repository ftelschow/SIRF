#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Load packages
suppressPackageStartupMessages(library(SIRF))
suppressPackageStartupMessages(library(SampleFields))
suppressPackageStartupMessages(library(RFT))
suppressPackageStartupMessages(library(tidyverse))

# Get the working directory
path_wd = args[6]

# Set the working path
setwd(path_wd)

# needed to make string without empty spaces
if(args[3] == "skewnessN"){
  transform = "skewness (normality)"
}else if(args[3] == "kurtosisN"){
  transform = "kurtosis (normality)"
}else{
  transform = args[3]
}

# Get the simulation file from source
source(paste(path_wd, "Article_simulation.R", sep = ""))

# Simulate the covering rate
Article_simulation( Model   = args[1], # "ModelA", "ModelB", "ModelC"
                    Msim    = as.numeric(args[2]), # number of simulations
                    transformation = transform, # Transformation
                    Nvec    = c( 50, 100, 200, 400, 800 ), # sample sizes considered
                    x       = seq( 0, 1, length.out = 150 ), # locations the process is evaluated at
                    level   = .95,  # level of simultaneous control
                    sim_tag = args[4],
                    path_wd = path_wd,
                    date    = args[5] )
