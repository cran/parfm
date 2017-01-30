################################################################################
#  Plot of objects of class 'select.parfm'                                     #
################################################################################
#                                                                              #
#                                                                              #
#                                                                              #
#   Date: December 22, 2011                                                    #
#   Last modification on: December  2, 2016                                    #
################################################################################

plot.select.parfm <- function(x, 
                              mar=c(2.5, 2, 1.5, .5),
                              ty = 'b',
                              ...){
  par(mfrow=c(1, 3))
  
  ### --- AIC --- ###
  par(mar=mar)
  plot(0,0, ty="n", xlab="", ylab="", main="AIC", xaxt="n",
    ylim=c(min(x$AIC, na.rm=TRUE) *  .9975,
           max(x$AIC, na.rm=TRUE) * 1.0025),
    xlim=c(.5, ncol(x$AIC) + .5), cex.lab=1.5)
  abline(v=1:ncol(x$AIC), col="grey")    
  
  mtext(c(none="No",
          gamma="Ga",
          ingau="IG",
          possta="PS",
          lognor="LN")[colnames(x$AIC)],
        side=1, at=1:ncol(x$AIC), padj=1)
  
  for (i in 1:nrow(x$AIC)) points(
    (1:ncol(x$AIC)), x$AIC[i, ],
    pch = 19 + i, cex = 1.5, ty = ty, bg = i)
  
  
  ### --- names --- ###
  par(mar=mar)
  plot(0:2, 0:2, xaxt = "n", yaxt = "n", bty = "n", ann = FALSE,
       ty = "n")
  
  legend("top", #c(.3, 1.7), c(1, 1.75),
         title = 'Baseline',
         c(exponential = "exponential",
           weibull = "Weibull", 
           inweibull = "inverse Weibull",
           gompertz = "Gompertz",
           loglogistic = "loglogistic", 
           lognormal = "logNormal",
           logskewnormal = "logSkewNormal")[rownames(x$AIC)],
         pch = {if(ty == 'l') NULL else 19 + 1:nrow(x$AIC)},
         pt.bg = 1:nrow(x$AIC),
         bg = "white", bty = "n", lty = ifelse(ty == 'p', 0, 1),
         ncol = 1, cex = 1.5, xjust = .5)
  
  legend("bottom", #c(0, 2), c(.25, 1), yjust=1,
         title = 'Frailty distribution',
         mapply(paste, 
                c(none="No",
                  gamma="Ga",
                  ingau="IG",
                  possta="PS",
                  lognor="LN")[colnames(x$AIC)],
                c(none="no frailty",
                  gamma="gamma",
                  ingau="inverse Gaussian",
                  possta="positive stable",
                  lognor="lognormal")[colnames(x$AIC)],
                sep=" = "),
         bg="white", bty="n",
         ncol=1, cex=1.5, xjust=.5)
  ### --- end names --- ###
  
  
 
  ### --- BIC --- ###
  par(mar=mar)
  plot(0,0, ty="n", xlab="", ylab="", main="BIC", xaxt="n",
    ylim=c(min(x$BIC, na.rm=TRUE) *  .9975,
           max(x$BIC, na.rm=TRUE) * 1.0025),
    xlim=c(.5, ncol(x$BIC) + .5), cex.lab=1.5)
  abline(v=1:ncol(x$BIC), col="grey")    
  
  mtext(c(none="No",
          gamma="Ga",
          ingau="IG",
          possta="PS",
          lognor="LN")[colnames(x$BIC)],
        side=1, at=1:ncol(x$BIC), padj=1)
  
  for (i in 1:nrow(x$BIC)) points(
    (1:ncol(x$BIC)), x$BIC[i, ],
    pch = 19 + i, cex = 1.5, ty = ty, bg = i)
}

