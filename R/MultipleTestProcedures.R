#------------------------------------------------------------------------------#
#                                                                              #
#     Functions to compute simultaneous confidence bands for functional data   #
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#  - scb_moments
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline included
# - add method description
#------------------------------------------------------------------------------#
#' This functions computes the true standard error depending on N for different
#' functionals of Gaussians
#'
#' @param alpha numeric in (0,1) the overall significance level to be achieved
#' @param pvals a vector of pvalues
#' @param correction string currently "bonferroni" or "holm" are options
#' @return Standard error under the assumption the data is Gaussian
#' @export
FWE_control <- function(alpha, pvals, correction){
  M = length(pvals)
  reject <- rep(FALSE, M)

  if(correction == "bonferroni"){
    reject = pvals < alpha / M
  }else if(correction == "holm"){
    Isort <- base::order(pvals, decreasing = FALSE)
    pvals_sort <- pvals[Isort]
    pvals_sort <- vapply( 1:M, function(p) (M + 1 - p) * pvals_sort[p], FUN.VALUE = 0.1 )
    reject[Isort[pvals_sort < alpha]] = TRUE
  }else if(correction == "sidak"){
    Isort <- base::order(pvals, decreasing = FALSE)
    pvals_sort <- pvals[Isort]
    alpha_m <- vapply( 1:M, function(p) 1-(1-alpha)^(1/(M+1-p)), FUN.VALUE = 0.1 )
    reject[Isort[pvals_sort < alpha_m]] = TRUE
  }else{
    stop("Method not implmenented. Choose 'bonferroni' for Bonferroni
         correction or 'BH' for Bonferroni-Holm Correction.")
  }
  return(reject)
}

#' @inheritParams FWE_control
#' @param correction string currently "BH" or "BY" are options
#' @return Standard error under the assumption the data is Gaussian
#' @export
FDR_control <- function(alpha, pvals, correction){
  # Get the number of the tests
  M = length(pvals)
  # Initialize the rejection vector
  reject <- rep(FALSE, M)
  # Save the order and sort the pvalues
  Isort      <- base::order(pvals, decreasing = FALSE)
  pvals.sort <- pvals[Isort]

  if(correction %in% c("BH", "BY")){
    if(correction == "BH"){
      C = rep(1, M)
    }else{
      C = cumsum(1/1:M)
    }
    # "correct" the pvalue
    pvals_sort <- vapply( 1:M, function(k) (M/k) * C[k] * pvals.sort[k], FUN.VALUE = 0.1 )
    # BH rule
    k_max      <- max(which(pvals_sort < alpha))
    # Reject all sorted k smaller than k_max
    reject[Isort[1:k_max]] = TRUE
  }else{
    stop("Method not implmenented. Choose 'bonferroni' for Bonferroni
         correction or 'BH' for Bonferroni-Holm Correction.")
  }
  return(reject)
}
