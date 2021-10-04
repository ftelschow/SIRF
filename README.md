# Statistical Inference for Random Fields

This repository contains a R-package implementing different methods for
simultaneous inference of random fields.

Especially, it contains the simulation code for different articles concerning
these questions.

## SIRF
 This package implements methods for Simultaneous Inference
 for Random Fields (SIRF).
 Currently, the main functionality is computation of simultaneous confidence
 bands for the population mean or moment based statistics of differentiable
 random field data over a rectangular domain in 1D or 2D.
 Options to estimate the quantile are a non-parametric bootstrap,
 a multiplier bootstrap or the GKF. The latter is the main novelty since
 it allows for fast computation of the quantiles even for higher dimensional
 domains.
