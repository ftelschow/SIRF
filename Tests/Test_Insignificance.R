  #-------------------------------------------------------------------------------
  # This file tests the preimage estimation
  #-------------------------------------------------------------------------------
  # Prepare workspace
  rm(list = ls())

  library(SampleFields)
  library(tidyverse)
  library(SIRF)
  library(mutoss)
  library(multtest)

  setwd("/home/fabian/Seafile/Code/Rpackages/SIRF/")

  source("Tests/Auxillary_fcns.R")
  #-------------------------------------------------------------------------------
  # Variables: General
  N     = 1e2
  alpha = 0.1
  betaN = 0.99
  k_fac = 2

  mcorrect = "sidak" # "holm" #

  # Variables: mean function
  mu_name = "1" # "3" #  "2" #
  p_mu2   = 0.05
  nx_true_low = 20
  nx_true_up  = 0
  nx_close    = 30
  nx_close2   = 30
  nx_far      = 0
  # variables: q estimation
  name       = "t" #"gauss" # "mboot" #
  SCoPEStype = "classical" # "selection" # "extraction" #
  mu1est     = "thickening" # NULL #
  truesigma  = FALSE

  kN = log(N) / k_fac


  mm = list(minus = rep(T, length(muvec)),#-nx_close2),
            plus = rep(T, length(muvec)))#-nx_close2))
  qest <- function(q) maxT_p(q, mm, df = N-1) - (1 - betaN)

  SCBquant = uniroot(qest, interval = c(-100, 100))
  kN = SCBquant$root
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
  if(SCoPEStype == "classical"){
    I0 = model$mu(model$x) == C
    I1 = model$mu(model$x) != C
  }

  # Get the SCoPES for the method
  scopes <- SCoPES(alpha = alpha, C = C, x = x, hatmu = hatmu,
                   hatsigma = hatsigma, tN = 1 / sqrt(N),
                   method = method,  R = R, mu = model$mu(model$x))

  plot_SCoPES(scopes = scopes, index_C = 1, mu = mu)

  # IVs = IVobs(scopes, method)
  # H0 = any(scopes$hatLC) | any(scopes$hatUC)
  #
  # IV = round(c(H0, IVs$IV0, IVs$IV1, IVs$IV1sim, IVs$IVk, IVs$IVkabs, IVs$IVobs), 3)
  # names(IV) <- c("H0", "IV0", "IV1", "IV1sim", "IVk", "IVkabs", "IVobs")
  # IVloc = round(IVs$IVloc, 3)

  if(SCoPEStype == "classical"){
    # Bonferroni-Holm-correction
    pvals <- apply(y$Y, 1, function(v) t.test(x = v,
                              alternative = "two.sided",
                              conf.level = 1-alpha)$p.value)
   # pvals_adjust = p.adjust(pvals, method = "sidak")
    detect_holm  = FWE_control(alpha, pvals, mcorrect)

    t_detect_h <- sum(detect_holm[I1])
    f_detect_h <- sum(detect_holm[I0])

    test = round(cbind(c(t_detect_h, f_detect_h), c(scopes$NtrueDetect, scopes$NfalseDetect)), 3)
    colnames(test) <- c(mcorrect, "SCoPES")
    rownames(test) <- c("true", "false")

    # IVloc = rbind(IVloc, pvals, pvals_adjust)
    # colnames(IVloc) <- x
    # rownames(IVloc) <- c("IVloc", "pvals", "pvals.adj")
  }else{
    test = round(c(scopes$NtrueDetect, scopes$NfalseDetect), 3)
    names(test) <- c("true", "false")
  }

  test
  # IV
  #IVloc
kN
