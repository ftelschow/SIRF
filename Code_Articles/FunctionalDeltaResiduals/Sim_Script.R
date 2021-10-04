#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

# Load packages
library(SCBfun)
library(SampleFields)
library(RFT)

# Get the working directory
path_wd = args[5]

# Set the working path
setwd(path_wd)

# Get the simulation file from source
source(paste(path_wd, "Article_simulation1.R", sep = ""))

# Simulate the covering rate
Article_simulation1( Model   = args[1], # "ModelA", "ModelB", "ModelC"
                     Msim    = as.numeric(args[2]), # number of simulations
                     Nvec    = c( 50, 100, 200, 400, 800 ), # sample sizes considered
                     x       = seq( 0, 1, length.out = 150 ), # locations the process is evaluated at
                     level   = .95,  # level of simultaneous control
                     sim_tag = args[3],
                     path_wd = path_wd,
                     date    = args[4] )
