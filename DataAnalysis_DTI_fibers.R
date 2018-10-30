###################################################################################
####
#### Final version of DTI fiber analysis, Aug 3 2018
####
###################################################################################
#################### Load required packages
library(SCBfda)
library(R.matlab)
library(gtools)
library(fda)
library(nortest)
require(fields)
rm(list=ls())

setwd("/media/sf_Linux/Research/Projects/2017_GKFinFDA")

########################## Load data ###########################################
data       <- readMat("Data/afq_fa.mat")
DTIdata    <- data$fa
fiberNames <- unlist(data$fgnames)
rm(data)

########################## Visualize data ######################################
level  = 0.95
Nbw    = 200
xeval  = seq(0,100, length.out=300)
x      = seq(0,100, length.out=100)
bw     = 100*seq( 0.015, 0.1, length.out=Nbw )

# pdfname <- paste( "Pics/", 100*level, "DTIdataVisual.pdf", sep="" )
# pdf(pdfname,  width=15, height=25)
# par(mfrow=c(7,4))
# par(mar=c(2.1,5.1,3.1,2.1) )
# for(k in 1:28){
#   #### Remove NA measurements
#   y1 <- NULL
#   y2 <- NULL
#
#   for(n in 1:15){
#     if(all(!is.na(DTIdata[, n,k]))){
#       y1 <- cbind(y1, DTIdata[, n,k])
#     }
#     if(all(!is.na(DTIdata[, n+15,k]))){
#       y2 <- cbind(y2, DTIdata[, n+15,k])
#     }
#   }
#
#   #### Plot the data with SCBs
#   plot( NULL,
#         xlim = c(min(x),max(x)),
#         ylim = c(0,1),
#         xlab = "location",
#         main = k,
#         ylab = fiberNames[k]
#   )
#   ## draw data
#   matlines( x, y1, lty=1, col="lightblue" )
#   matlines( x, y2, lty=1, col="pink" )
#   legend("topleft",
#          legend = c("control","patient"),
#          lty=c(1,1),
#          col= c("lightblue", "pink"), bty="n",
#          cex=1.3,
#          yjust=1,
#          xjust=0.5
#   )
#   ## calculate SCB & draw
#   SCB.control <- scb_mean( y1, level=level, method="tGKF", param_method=NULL )
#   SCB.patient <- scb_mean( y2, level=level, method="tGKF", param_method=NULL )
#   SCB2.control <- scb_mean( y1, level=level, method="Bootstrapt", param_method=NULL )
#   SCB2.patient <- scb_mean( y2, level=level, method="Bootstrapt", param_method=NULL )
#
#   matlines( x, SCB.control$hatmean, col="blue" )
#   matlines( x, cbind(SCB.control$scb$lo, SCB.control$scb$up), col="blue", lty=2 )
#   matlines( x, cbind(SCB2.control$scb$lo, SCB2.control$scb$up), col="blue", lty=2 )
#   matlines( x, SCB.patient$hatmean, col="red" )
#   matlines( x, cbind(SCB.patient$scb$lo, SCB.patient$scb$up), col="red", lty=2 )
#   matlines( x, cbind(SCB2.patient$scb$lo, SCB2.patient$scb$up), col="orange", lty=2 )
# }
# dev.off()
# ########################### Q-Q-Plot ###########################
# AD.pvals <- array(NA, dim=c(28,10,3))
#
# pdfname <- paste( "Pics/", 100*level, "DTIdataQQ.pdf", sep="" )
# pdf(pdfname,  width=15, height=25)
# par(mfrow=c(10,3))
# for(k in 1:28){
#   #### Remove NA measurements
#   y1 <- NULL
#   y2 <- NULL
#
#   for(n in 1:15){
#     if(all(!is.na(DTIdata[, n,k]))){
#       y1 <- cbind(y1, DTIdata[, n,k])
#     }
#     if(all(!is.na(DTIdata[, n+15,k]))){
#       y2 <- cbind(y2, DTIdata[, n+15,k])
#     }
#   }
#
#   res1 = y1 - rowMeans(y1)
#   res2 = y2 - rowMeans(y2)
#   subTimes <- seq(1,100, length.out=10)
#   for(i in 1:10){
#     nt = subTimes[i]
#     y = res1[nt,]
#     qqnorm(y)
#     qqline(y, col = "red")
#     AD.pvals[k,i,1] <- ad.test(y)$p
#     y = res2[nt,]
#     qqnorm(y)
#     qqline(y, col = "red")
#     AD.pvals[k,i,2] <- ad.test(y)$p
#     y = c(res1[nt,], res2[nt,])
#     qqnorm(y)
#     qqline(y, col = "red")
#     AD.pvals[k,i,3] <- ad.test(y)$p
#   }
# }
# dev.off()
#
# thresh <- 0.95
# length(which(AD.pvals[,,1]<thresh))/length(AD.pvals[,,1])
# length(which(AD.pvals[,,2]<thresh))/length(AD.pvals[,,2])
# length(which(AD.pvals[,,3]<thresh))/length(AD.pvals[,,3])
#
#
# #####################################################################
# #####################################################################
# #############################  Testing  #############################
# #####################################################################
# #####################################################################
#
# ##### Compute weights for local linear smoother
# SmoothWeights = array( NA, dim=c(length(xeval), length(x), Nbw) )
# for( j in 1:Nbw ){
#   #### Compute weights for local linear regression
#   SmoothWeights[,,j]   <- as.matrix( locpol::locLinWeightsC( x = x, xeval = xeval, bw = bw[j], kernel = locpol::gaussK )$locWeig )
# }
#
#
# pdfname <- paste( "Pics/", 100*level,"bwLim", 10*max(bw), "DTIresults.pdf", sep="" )
# pdf(pdfname,  width=15, height=25)
# par(mfrow=c(7,4))
# for(k in 1:28){
#   #### Remove NA measurements
#   y1 <- NULL
#   y2 <- NULL
#
#   for(n in 1:15){
#     if(all(!is.na(DTIdata[, n,k]))){
#       y1 <- cbind(y1, DTIdata[, n,k])
#     }
#     if(all(!is.na(DTIdata[, n+15,k]))){
#       y2 <- cbind(y2, DTIdata[, n+15,k])
#     }
#   }
#   ## Transform the samples into a scale field
#   sm.y1 = scaleField(Y=y1, x=x, xeval=xeval, h=bw, Weights = SmoothWeights)
#   sm.y2 = scaleField(Y=y2, x=x, xeval=xeval, h=bw, Weights = SmoothWeights)
#   scalebds <- scb_meandiff( Y1 = sm.y1, Y2 = sm.y2, level = level, method="tGKF" )
#
#   # test
#   Test.results <- 1*( !( (scalebds$scb$lo<0) & (scalebds$scb$up>0) ) )
#   # Plotting the confidence bands
#   fields::image.plot(Test.results, xlab = "location", ylab = "bandwidth",
#                    #  xlim = range(x),
#                   #   ylim = range( bw ), zlim=c(-0.001,1.5),
#                      main = paste("DTIfiber ", k, " ScaleBands Test ",100*level,"%", sep=""))
#   legend("topleft",
#          legend = c("Accept", "Reject"),
#          lty=c(2,2),
#          col= c("blue", "red"), bty="n",
#          cex=1.3,
#          yjust=1,
#          xjust=0.5
#   )
#   print(k)
# }
# dev.off()
#
# ############### No scale Space
# pdfname <- paste( "Pics/", 100*level,"DTIresultsNoScale.pdf", sep="" )
# pdf(pdfname,  width=15, height=25)
# par(mfrow=c(7,4))
# for(k in 1:28){
#   #### Remove NA measurements
#   y1 <- NULL
#   y2 <- NULL
#
#   for(n in 1:15){
#     if(all(!is.na(DTIdata[, n,k]))){
#       y1 <- cbind(y1, DTIdata[, n,k])
#     }
#     if(all(!is.na(DTIdata[, n+15,k]))){
#       y2 <- cbind(y2, DTIdata[, n+15,k])
#     }
#   }
#   nSamp1 = dim(y1)[2]
#   nSamp2 = dim(y2)[2]
#   scb <- SCBglm( x = times, y = cbind(y1, y2), X=cbind( c(rep(1,nSamp1), rep(0,nSamp2)), c(rep(0,nSamp1), rep(1,nSamp2))), c=c(1,-1), xlim = c(0,Tmax), bw=0.01, level = level )
#
#   plot( NULL,
#         xlim=range(x), ylim=range(scb$scb),
#         xlab = "location", ylab = fiberNames[k],
#         main = k
#   )
#   matlines( seq(0,1,length.out=400), scb$scb, lty=2, col=2 )
#   lines( seq(0,1,length.out=400), scb$hatmean, lwd=1.5, col=2)
#   abline(h=0)
#   print(k)
# }
# dev.off()
#
# #### Random Pairings
# ############### No scale Space
# pdfname <- paste( "Pics/", 100*level,"DTIresultsRandomPairsNoScale.pdf", sep="" )
# pdf(pdfname,  width=15, height=25)
# par(mfrow=c(7,4))
# for(k in 1:28){
#   #### Remove NA measurements
#   y1 <- NULL
#   y2 <- NULL
#
#   for(n in 1:15){
#     if(all(!is.na(DTIdata[, n,k]))){
#       y1 <- cbind(y1, DTIdata[, n,k])
#     }
#     if(all(!is.na(DTIdata[, n+15,k]))){
#       y2 <- cbind(y2, DTIdata[, n+15,k])
#     }
#   }
#
#   nSamp1 = dim(y1)[2]
#   nSamp2 = dim(y2)[2]
#   if( all( c(nSamp1,nSamp2)==c(15,15) ) ){
#     print(c(k,nSamp1,nSamp2))
#     # scb <- SCBglm( x = times, y = cbind(y1, y2), X=cbind( c(rep(1,nSamp1), rep(0,nSamp2)), c(rep(0,nSamp1), rep(1,nSamp2))), c=c(1,-1), xlim = c(0,Tmax), bw=0.01, level = level )
#     #
#     # plot( NULL,
#     #       xlim=range(times), ylim=range(scb$scb),
#     #       xlab = "location", ylab = fiberNames[k],
#     #       main = k
#     # )
#     # matlines( seq(0,1,length.out=400), scb$scb, lty=2, col=2 )
#     # lines( seq(0,1,length.out=400), scb$hatmean, lwd=1.5, col=2)
#     # abline(h=0)
#     # print(k)
#   }
# }
# dev.off()


#################################################################################
#################################################################################
#############################  Plots for the paper  #############################
#################################################################################
#################################################################################

##### Compute weights for local linear smoother
SmoothWeights = array( NA, dim=c(length(xeval), length(x), Nbw) )
for( j in 1:Nbw ){
  #### Compute weights for local linear regression
  SmoothWeights[,,j]   <- as.matrix( locpol::locLinWeightsC( x = x, xeval = xeval, bw = bw[j], kernel = locpol::gaussK )$locWeig )
}


pdfname <- paste( "Pics/", 100*level,"bwLim", 10*max(bw), "DTIresults.pdf", sep="" )
pdf(pdfname,  width=0.7*13, height=0.7*9)
par(mfrow=c(3,3))
par(mar=c(2.1,5.1,3.1,2.1) )

  #### Remove NA measurements
  y2.c <- y8.c <- y15.c <- NULL
  y2.p <- y8.p <- y15.p <- NULL

  for(n in 1:15){
    if(all(!is.na(DTIdata[, n,2]))){
      y2.c <- cbind(y2.c, DTIdata[, n,2])
    }
    if(all(!is.na(DTIdata[, n+15,2]))){
      y2.p <- cbind(y2.p, DTIdata[, n+15,2])
    }
  }
  for(n in 1:15){
    if(all(!is.na(DTIdata[, n,8]))){
      y8.c <- cbind(y8.c, DTIdata[, n,8])
    }
    if(all(!is.na(DTIdata[, n+15,2]))){
      y8.p <- cbind(y8.p, DTIdata[, n+15,8])
    }
  }
  for(n in 1:15){
    if(all(!is.na(DTIdata[, n,15]))){
      y15.c <- cbind(y15.c, DTIdata[, n,15])
    }
    if(all(!is.na(DTIdata[, n+15,2]))){
      y15.p <- cbind(y15.p, DTIdata[, n+15,15])
    }
  }

  ## Plot the data and its mean
  # Fiber 2
  plot( NULL,
        xlim = c(min(x),max(x)),
        ylim = c(0.2,0.8),
        xlab = "Location",
        ylab = fiberNames[2],
        cex.lab = 1.3,
        cex.axis= 1.3
  )
  ## draw data
  matlines( x, y2.c, lty=1, col="lightblue" )
  matlines( x, y2.p, lty=1, col="pink" )
  lines(x, rowMeans(y2.c), lty=1, lwd=1.5, col="blue")
  lines(x, rowMeans(y2.p), lty=1, lwd=1.5, col="red")

  # Fiber 8
  plot( NULL,
        xlim = c(min(x),max(x)),
        ylim = c(0.2,0.8),
        xlab = "Location",
        ylab = fiberNames[8],
        cex.lab = 1.3,
        cex.axis= 1.3
  )
  ## draw data
  matlines( x, y8.c, lty=1, col="lightblue" )
  matlines( x, y8.p, lty=1, col="pink" )
  lines(x, rowMeans(y8.c), lty=1, lwd=1.5, col="blue")
  lines(x, rowMeans(y8.p), lty=1, lwd=1.5, col="red")
  legend("topleft",
         legend = c("control","patient"),
         lty=c(1,1),
         col= c("lightblue", "pink"), bty="n",
         cex=1.3,
         yjust=1,
         xjust=0.5
  )
  # Fiber 15
  plot( NULL,
        xlim = c(min(x),max(x)),
        ylim = c(0.2,0.8),
        xlab = "Location",
        ylab = fiberNames[15],
        cex.lab = 1.3,
        cex.axis= 1.3
  )
  ## draw data
  matlines( x, y15.c, lty=1, col="lightblue" )
  matlines( x, y15.p, lty=1, col="pink" )
  lines(x, rowMeans(y15.c), lty=1, lwd=1.5, col="blue")
  lines(x, rowMeans(y15.p), lty=1, lwd=1.5, col="red")

  ## Plot and compute confidence bands of the difference
  # fiber 2
  scb_tGKF   = scb_meandiff( Y1=y2.c, Y2=y2.p, level = level, method="tGKF" )
  scb_mboott = scb_meandiff( Y1=y2.c, Y2=y2.p, level = level, method="MultiplierBootstrapt" )
  plot( NULL,
        xlim = c(min(x),max(x)),
        ylim =  c(-0.2,0.2),
        xlab = "Location",
        ylab = fiberNames[2],
        cex.lab = 1.3,
        cex.axis= 1.3
  )
  abline(h = 0, col=1)
  lines( x, scb_tGKF$hatmean, lty=1, lwd=1.5, col="grey3")
  lines( x, scb_mboott$scb$lo, lty=2, lwd=1.5, col="blue" )
  lines( x, scb_mboott$scb$up, lty=2, lwd=1.5, col="blue" )
  lines( x, scb_tGKF$scb$lo, lty=2, lwd=1.5, col="red" )
  lines( x, scb_tGKF$scb$up, lty=2, lwd=1.5, col="red" )


  # fiber 8
  scb_tGKF   = scb_meandiff( Y1=y8.c, Y2=y8.p, level = level, method="tGKF" )
  scb_mboott = scb_meandiff( Y1=y8.c, Y2=y8.p, level = level, method="MultiplierBootstrapt" )
  plot( NULL,
        xlim = c(min(x),max(x)),
        ylim =  c(-0.2,0.2),
        xlab = "Location",
        ylab = fiberNames[8],
        cex.lab = 1.3,
        cex.axis= 1.3
  )
  abline(h = 0, col=1)
  lines( x, scb_tGKF$hatmean, lty=1, lwd=1.5, col="grey3")
  lines( x, scb_mboott$scb$lo, lty=2, lwd=1.5, col="blue" )
  lines( x, scb_mboott$scb$up, lty=2, lwd=1.5, col="blue" )
  lines( x, scb_tGKF$scb$lo, lty=2, lwd=1.5, col="red" )
  lines( x, scb_tGKF$scb$up, lty=2, lwd=1.5, col="red" )

  legend("topleft",
         legend = c("tGKF","Multiplier-t"),
         lty=c(1,1),
         col= c("red", "blue"), bty="n",
         cex=1.3,
         yjust=1,
         xjust=0.5
  )
  # fiber 15
  scb_tGKF   = scb_meandiff( Y1=y15.c, Y2=y15.p, level = level, method="tGKF" )
  scb_mboott = scb_meandiff( Y1=y15.c, Y2=y15.p, level = level, method="MultiplierBootstrapt" )
  plot( NULL,
        xlim = c(min(x),max(x)),
        ylim = c(-0.2,0.2),
        xlab = "Location",
        ylab = fiberNames[15],
        cex.lab = 1.3,
        cex.axis= 1.3
  )
  abline(h = 0, col=1)
  lines( x, scb_tGKF$hatmean, lty=1, lwd=1.5, col="grey3")
  lines( x, scb_mboott$scb$lo, lty=2, lwd=1.5, col="blue" )
  lines( x, scb_mboott$scb$up, lty=2, lwd=1.5, col="blue" )
  lines( x, scb_tGKF$scb$lo, lty=2, lwd=1.5, col="red" )
  lines( x, scb_tGKF$scb$up, lty=2, lwd=1.5, col="red" )


  #### Plot and compute scale process confidence bands
  ## Fiber 2
  sm.y1 = scaleField(Y=y2.c, x=x, xeval=xeval, h=bw, Weights = SmoothWeights)
  sm.y2 = scaleField(Y=y2.p, x=x, xeval=xeval, h=bw, Weights = SmoothWeights)
  scalebds <- scb_meandiff( Y1 = sm.y1, Y2 = sm.y2, level = level, method="tGKF" )

  # test
  Test.results <- 1*( !( (scalebds$scb$lo<0) & (scalebds$scb$up>0) ) )
  # Plotting the confidence bands
  image(Test.results, xlab = "Location", ylab = "Bandwidth", axes=F, col=terrain.colors(100), cex.lab=1.3 )
  axis(1,  at=seq(0,1,length.out=5), labels=seq(0,100, length.out=5), cex.axis=1.3 )
  axis(2,  at=seq(0,1,length.out=6), labels=seq(1.5,10, length.out=6), cex.axis=1.3 )
  legend("topleft",
         legend = c("Accept", "Reject"),
         lty=c(2,2), lwd = c(2,2),
         col= c("green", "white"), bty="n",
         cex=1.3,
         yjust=1,
         xjust=0.5
  )
  ## fiber 8
  sm.y1 = scaleField(Y=y8.c, x=x, xeval=xeval, h=bw, Weights = SmoothWeights)
  sm.y2 = scaleField(Y=y8.p, x=x, xeval=xeval, h=bw, Weights = SmoothWeights)
  scalebds <- scb_meandiff( Y1 = sm.y1, Y2 = sm.y2, level = level, method="tGKF" )

  # test
  Test.results <- 1*( !( (scalebds$scb$lo<0) & (scalebds$scb$up>0) ) )
  # Plotting the confidence bands
  image(Test.results, xlab = "Location", ylab = "Bandwidth", axes=F, col=terrain.colors(100), cex.lab=1.3 )
  axis(1,  at=seq(0,1,length.out=5), labels=seq(0,100, length.out=5), cex.axis=1.3 )
  axis(2,  at=seq(0,1,length.out=6), labels=seq(1.5,10, length.out=6), cex.axis=1.3 )
  legend("topleft",
         legend = c("Accept", "Reject"),
         lty=c(2,2), lwd = c(2,2),
         col= c("green", "white"), bty="n",
         cex=1.3,
         yjust=1,
         xjust=0.5
  )
  ## fiber 15
  sm.y1 = scaleField(Y=y15.c, x=x, xeval=xeval, h=bw, Weights = SmoothWeights)
  sm.y2 = scaleField(Y=y15.p, x=x, xeval=xeval, h=bw, Weights = SmoothWeights)
  scalebds <- scb_meandiff( Y1 = sm.y1, Y2 = sm.y2, level = level, method="tGKF" )

  # test
  Test.results <- 1*( !( (scalebds$scb$lo<0) & (scalebds$scb$up>0) ) )
  # Plotting the confidence bands
  image(Test.results, xlab = "Location", ylab = "Bandwidth", axes=F, col=terrain.colors(100), cex.lab=1.3 )
  axis(1,  at=seq(0,1,length.out=5), labels=seq(0,100, length.out=5), cex.axis=1.3 )
  axis(2,  at=seq(0,1,length.out=6), labels=seq(1.5,10, length.out=6), cex.axis=1.3 )
  legend("topleft",
         legend = c("Accept", "Reject"),
         lty=c(2,2), lwd = c(2,2),
         col= c("green", "white"), bty="n",
         cex=1.3,
         yjust=1,
         xjust=0.5
  )
dev.off()
