#------------------------------------------------------------------------------#
#                                                                              #
#     Miscalaneous functions                                                   #
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#  - scaleField
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline NOT included
#------------------------------------------------------------------------------#
#' This functions computes the standard
#'
#' @param x The considered transformation of the moments
#' @param scb Vector grid for evaluation of the smoothed field
#' @param scb2 Vector grid for evaluation of the smoothed field
#' @param true Vector grid for evaluation of the smoothed field
#' @param title Vector grid for evaluation of the smoothed field
#' @return Standard error under the assumption the data is Gaussian
#' @export
plot_scb <- function(x, scb, scb2 = NULL, true = NULL, s1name = "est Sample1",
                     s2name = "est Sample2", title = "", legend.position = "none"){
  N1 = dim(scb$Y)[2]

  if(is.null(true)){
    data = cbind(x, scb$scb)
    data = as_tibble(data)
  }else{
    data = cbind(x, scb$scb, true)
    colnames(data)[5] = "true"
    data = as_tibble(data)
  }

  if(!is.null(scb2)){
    data2 = cbind(x, scb2$scb)
    data2 = as_tibble(data2)
    N2 = dim(scb2$Y)[2]
  }

  sLegend <- 18
  sText   <- 18
  Sylab   <- 1.5
  Sxlab   <- 1.5
  sLine   <- 1.5
  sPch    <- 4

  theme2 <- theme(legend.position = legend.position,
                  legend.title = element_blank(),
                  legend.text  = element_text(color = "black", size = sText, face = "plain"),
                  legend.direction = "horizontal",
                  legend.key = element_rect(colour = "transparent", fill = "transparent"),
                  legend.background = element_blank(),
                  plot.title   = element_text(face = "bold", size = sText),
                  axis.text.x  = element_text(color = "black", size = sText, face = "plain"),
                  axis.text.y  = element_text(color = "black", size = sText, face = "plain"),
                  axis.title.x = element_text(color = "black", size = sText, face = "plain"),
                  axis.title.y = element_text(color = "black", size = sText, face = "plain"))

  if(is.null(true) & is.null(scb2)){
    ggplot(data) +
      geom_line(aes(x = x, y = est, col = s1name), size = 1.2) +
      geom_ribbon(aes(x = x, ymin = SCB.low, ymax = SCB.up, fill = s1name), alpha = 0.25) +
      xlab( "" ) + ylab( "" ) +
      ggtitle(title) +
      theme2
  }else if(!is.null(true) & is.null(scb2)){
    ggplot(data) +
      geom_line(aes(x = x, y = est, col = s1name), size = 1.2) +
      geom_line(aes(x = x, y = true), size = 1.2) +
      geom_ribbon(aes(x = x, ymin = SCB.low, ymax = SCB.up, fill = s1name), alpha = 0.25) +
      xlab( "" ) + ylab( "" ) +
      ggtitle(title) +
      theme2
  }else if(is.null(true) & !is.null(scb2)){
    ggplot(data) +
      geom_line(aes(x = x, y = est, col = s1name), size = 1.2) +
      geom_ribbon(aes(x = x, ymin = SCB.low, ymax = SCB.up, fill = s1name), alpha = 0.25) +
      geom_line(aes(x = x, y = est, col = s2name), size = 1.2, data = data2) +
      geom_ribbon(aes(x = x, ymin = SCB.low, ymax = SCB.up, fill = s2name), data = data2, alpha = 0.25) +
      xlab( "" ) + ylab( "" ) +
      ggtitle(title) +
      theme2
  }else if(!is.null(true) & !is.null(scb2)){
    ggplot(data) +
      geom_line(aes(x = x, y = true), size = 1.2) +
      geom_line(aes(x = x, y = est, col = s1name), size = 1.2) +
      geom_ribbon(aes(x = x, ymin = SCB.low, ymax = SCB.up, fill = s1name), alpha = 0.25) +
      geom_line(aes(x = x, y = est, col = s2name), size = 1.2, data = data2) +
      geom_ribbon(aes(x = x, ymin = SCB.low, ymax = SCB.up, fill = s2name), data = data2, alpha = 0.25) +
      xlab( "" ) + ylab( "" ) +
      ggtitle(title) +
      theme2
  }


}

plot_covSim <- function(covSim, title = "", legend.position = "none"){
  lvl      <- covSim$level
  covRates <- covSim$rates
  Msim     <- covSim$Msim

  target <- sqrt( lvl * ( 1 - lvl ) / Msim ) * c(-qnorm(lvl/2 + 1/2), 0, qnorm(lvl/2 + 1/2)) + lvl
  xLab <- rownames(covRates)
  yLab <- colnames(covRates)
  covs <- as_tibble(covRates, rownames = "N") %>%
          reshape2::melt(., id.vars = "N", variable = "Method", value.name = "CovRate") %>%
          as_tibble() %>% mutate_if(is.character, as.numeric)

  # Plot the Covering Rates by Method
  print( ggplot(covs, aes(N, CovRate, group = Method, col = Method)) +
         geom_point() + geom_line() +
         xlab( "Sample Size [N]" ) + ylab( "Covering Rate" ) +
         ggtitle( title ) +
         geom_hline( yintercept = target[2] ) +
         geom_hline( yintercept = target[1], linetype = "dashed" ) +
         geom_hline( yintercept = target[3], linetype = "dashed" ) )
}
