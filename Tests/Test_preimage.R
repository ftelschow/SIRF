#-------------------------------------------------------------------------------
# This file tests the preimage estimation
#-------------------------------------------------------------------------------
# Prepare workspace
rm(list = ls())

library(SampleFields)
library(tidyverse)
library(SIRF)

source("Tests/Auxillary_fcns.R")
#-------------------------------------------------------------------------------
# Variables: General
N     = 20
alpha = 0.1

# Variables: mean function
muu   = c(0, -0.1, 0.1, 0.3)
muvec = c(rep(muu[1], 5), rep(muu[2], 5), rep(muu[3], 4), rep(muu[4], 6))
#muvec = c(1, rep(0, 6))

# variables: q estimation
name       = "t" #  "gauss" # "gauss" # "mboot" #
SCoPEStype = "extraction" #  "classical" #   "selection" #
mu1est     =  "thickening" # NULL #
kN         = log(N) / 2
truesigma  = TRUE

# Variables for the tubes
B = c(0, 3)

#-------------------------------------------------------------------------------
# Generate data
y = generateData(N, muvec, B, truesigma, SCoPEStype)
x = y$x; hatmu = y$hatmu; hatsigma = y$hatsigma; R = y$R; model = y$model;
C = y$C; tN = y$tN
#-------------------------------------------------------------------------------
# Generate method list
method = method_gen(name, SCoPEStype, mu1est, N, kN, R )


#-------------------------------------------------------------------------------
# Main test area
hatmu1Ctrue = PreimageC(C = C, hatmu = model$mu(model$x), hatsigma = hatsigma,
                        tN = tN, kN = 0, method = method$SCoPEStype)

hatmu1Cest = PreimageC(C = C, hatmu = hatmu, hatsigma = hatsigma,
                       tN = tN, kN = kN, method = method$SCoPEStype)

test <- c(sum(hatmu1Ctrue$minus), sum(hatmu1Cest$minus))
names(test) <- c("true", "est")
test
