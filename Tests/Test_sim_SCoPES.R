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
# Variables General
Msim  = 1e3
# Variables: General
N     = 5e1
alpha = 0.1

kNtype  = "SCB"
betaN = 0.1
k_fac = 2

# Variables: mean function
mu_name = "1" # "3" #  "2" #
p_mu2   = 0.05
nx_true_low = 50
nx_true_up  = 5
nx_close    = 25
nx_close2   = 35
nx_far      = 5
# variables: q estimation
name       = "t" # "gauss" # "mboot" #
SCoPEStype = "classical" # "extraction" #     "selection" #
mu1est     = "thickening" # NULL #
truesigma  = FALSE



# Variables for the tubes
Delta = 0.3

if(mu_name == "1"){
  B = c(0, 3)
  muu   = c(0, -0.3, 0.2, 3.5, 3)
  muvec = c(rep(muu[1], nx_true_low), rep(muu[2], nx_close2), rep(muu[3], nx_close), rep(muu[4], nx_far), rep(muu[5], nx_true_up))
  muvec = muvec[sample(1:length(muvec), replace = FALSE)]
}else if(mu_name == "2"){
  B = c(0, 3)
  muvec = c(p_mu2, rep(0, 20))
}else{
  muu   = c(0, B[1]-0.1, B[1]+0.1, B[2]-0.1, B[2]+0.1)
  muvec = c(rep(muu[1], 3000), rep(muu[2], 200), rep(muu[3], 200), rep(muu[4], 200), rep(muu[5], 200))
  B = c(-Delta, Delta)
}
kNold = log(N) / k_fac

mm = list(minus = rep(T, length(muvec)), plus = rep(T, length(muvec)))
qest <- function(q) maxT_p(q, mm, df = N-1) - 1 + betaN
SCBquant = uniroot(qest, interval = c(-50, 50))

c(kNold, SCBquant$root)
if(kNtype == "SCB"){
  kN = SCBquant$root
}else{
  kN = kNold
}

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
test = sim_SCoPES(Msim = Msim, N = N, alpha = alpha, C = C,
                  q_method = method, model = model)

results <- c(test$coverage)
names(results) <- c("Coverage")
results

if(SCoPEStype == "extraction"){
  index_minus = model$mu(model$x) == C[, 1]
  index_plus = model$mu(model$x) == C[, 2]

  par(mfrow = c(1,2))
  hist(apply(test$mu1C$minus[index_minus,], 2, sum),
       xlab = "Number Detected C-",
       main = paste("True = ", sum(index_minus), sep = "") )

  hist(apply(test$mu1C$plus[index_plus,], 2, sum),
       xlab = "Number Detected C+",
       main = paste( "True = ", sum(index_plus), sep = "") )

}else if(SCoPEStype == "classical"){
  index_minus = model$mu(model$x) == C
  par(mfrow = c(1,1))
  hist(apply(test$mu1C$minus[index_minus,], 2, sum),
       xlab = "Number Detected C",
       main = paste("True = ", sum(index_minus), sep = "") )
}

par(mfrow = c(1,1))

results <- c(test$coverage)
names(results) <- c("Coverage")
results

c(kNold, kN)

rowMeans(test$NDetect)
rowMeans(test$NDetect_holm)
rowMeans(test$NDetect_sidak)
