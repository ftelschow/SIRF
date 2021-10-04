#-------------------------------------------------------------------------------
#
#       Functions producing nice plots for confidence bands result
#
#-------------------------------------------------------------------------------
require(tidyverse)
require(reshape2)

plot_scb_cov <- function( covRates, Msim, lvl = 0.95, title = "Simultaneous Confidence Bands", ... ){
  target <- sqrt( lvl * ( 1 - lvl ) / Msim ) * c(-qnorm(lvl/2 + 1/2), 0, qnorm(lvl/2 + 1/2)) + lvl
  xLab <- rownames(covRates)
  yLab <- colnames(covRates)
  covs <- as_tibble(covRates, rownames = "N") %>%
                  melt(., id.vars = "N", variable = "Method", value.name = "CovRate") %>%
                              as_tibble() %>% mutate_if(is.character, as.numeric)
  
  # Plot the Covering Rates by Method
  tmp = ggplot(covs, aes(N, CovRate, group = Method, col = Method)) +
    geom_point() + geom_line() +
    xlab( "Sample Size [N]" ) + ylab( "Covering Rate" ) +
    ggtitle( title ) +
    geom_hline( yintercept = target[2] ) +
    geom_hline( yintercept = target[1], linetype = "dashed" ) +
    geom_hline( yintercept = target[3], linetype = "dashed" )
  return(tmp)
}

plot_scb <- function( scb,
                      title = "Simultaneous Confidence Bands",
                      xlab = "location",
                      ylab = "value",... ){

  require(reshape2)
  require(tidyverse)

  #locs   = scb$locations
  #values = scb$values
  values = scb$scb
  rownames(values) <- Y1$locations


  SCB <- as_tibble(values, rownames = "x") %>%
          melt(., id.vars = "x", variable = "Legend", value.name = "value") %>%
              as_tibble()

  # Plot the Covering Rates by Method
  ggplot(SCB, aes(x, value, group = Legend, col = Legend, linetype = Legend)) + geom_line() +
    scale_linetype_manual(values = c(2, 1, 2) ) +
    scale_color_manual(values = c('#9999CC', '#CC6666', '#9999CC')) +
    xlab( xlab ) + ylab( ylab ) + ggtitle( title )
}
