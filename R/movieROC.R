movieROC <- function(x, ...) {
  UseMethod("movieROC")
}

movieROC.groc <- function(x, fpr = NULL, h = 1, histogram = FALSE, breaks = 15, reduce = TRUE, completeROC = FALSE, videobar = TRUE, file = "animation1.gif", save = TRUE, legends = FALSE, speedcorrection = FALSE, tpause = 1, interval = 0.2, ani.width = ifelse(reduce, 830, 500), ani.height = ifelse(reduce, 300, 750), cex.lab = 2.5, cex.axis = 1.75, cex.main = 2.25+as.numeric(reduce), xlim = NULL, ylim = NULL, cex.point = 1.5, lwd.curve = 2, mar = NULL, lim.density = 0.01, xlab = "Marker",  main.density = "Density functions", col.controlscases = c('green4','red'), col.curve = 'black', col.threshold = 'blue', ...){

  histogram <- histogram; reduce <- reduce; breaks <- breaks; h <- h; completeROC <- completeROC; legends <- legends; reduce <- reduce; videobar <- videobar; tpause <- tpause
  xlim <- xlim; ylim <- ylim; cex.point <- cex.point; lwd.curve <- lwd.curve; cex.lab <- cex.lab; cex.axis <- cex.axis; cex.main <- cex.main; mar <- mar; lim.density <- lim.density; xlim <- xlim; xlab <- xlab;  main.density <- main.density; col.controlscases <- col.controlscases; col.curve <- col.curve; col.threshold <- col.threshold

  movie <- function(x, fpr){

    t <- x$t;

    if(is.null(fpr)){
      if(length(t) < 150 ) fpr <- t else fpr <- seq(min(t), max(t), length.out=100)
    }

    B <- length(fpr)

    if(videobar==TRUE){
      cat("\nProgress bar: Construction of GIF with ", B, " thresholds. \n", sep = "")
      bar <- txtProgressBar(min = 0, max = B, style = 3)
    }

    if(speedcorrection){

      output.predict0 <- predict(x, FPR = 0)$ClassSubsets
      output.predict <- predict(x, FPR = fpr[2])$ClassSubsets
      output.predict1 <- predict(x, FPR = 1)$ClassSubsets
      X <- sort(c(x$controls, x$cases))
      if(x$side == "right"){
        prop <- sum(X > output.predict[1] & X < output.predict0[1])/length(X)
        C.speed <- seq(output.predict0[1], output.predict[1], length.out = prop*length(fpr))
      }else if(x$side == "left"){
        prop <- sum(X > output.predict0[2] & X < output.predict[2])/length(X)
        C.speed <- seq(output.predict0[2], output.predict[2], length.out = prop*length(fpr))
      }
      if(x$side == "both"){
        prop <- sum((X > output.predict0[1,2] & X < output.predict[1,2])| (X > output.predict[2,1] & X < output.predict0[2,1]))/length(X)
        XL.speed <- seq(output.predict0[1,2], output.predict[1,2], length.out = prop*length(fpr))
        XU.speed <- seq(output.predict0[2,1], output.predict[2,1], length.out = prop*length(fpr))
      }else{
        prop <- sum(X > output.predict[1] & X < output.predict[2])/length(X)
        XL.speed <- seq(output.predict0[1], output.predict[1], length.out = prop*length(fpr))
        XU.speed <- seq(output.predict0[2], output.predict[2], length.out = prop*length(fpr))
      }

      sapply(c(1, seq(1+length(fpr),length(fpr)+prop*length(fpr),1), 2:length(fpr)), function(i){
        if(i <= length(fpr)){
          FPR <- fpr[i]
          plot.buildROC(x, FPR=FPR, h=h, build.process=TRUE, histogram=histogram, breaks=breaks, completeROC=completeROC, legends=legends, reduce=reduce, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, xlim = xlim, ylim = ylim, cex.point = cex.point, lwd.curve = lwd.curve, mar = mar, lim.density = lim.density, xlab=xlab, main.density = main.density, col.controlscases = col.controlscases, col.curve = col.curve, col.threshold = col.threshold)
        }else{
          if(x$side %in% c("right","left")){
            plot.buildROC(x, FPR=NULL, C=C.speed[i-length(fpr)], h=h, build.process=TRUE, histogram=histogram, breaks=breaks, completeROC=completeROC, legends=legends, reduce=reduce, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, xlim = xlim, ylim = ylim, cex.point = cex.point, lwd.curve = lwd.curve, mar = mar, lim.density = lim.density, xlab=xlab, main.density = main.density, col.controlscases = col.controlscases, col.curve = col.curve, col.threshold = col.threshold)
          }else{
            plot.buildROC(x, FPR=NULL, XL=XL.speed[i-length(fpr)], XU=XU.speed[i-length(fpr)], h=h, build.process=TRUE, histogram=histogram, breaks=breaks, completeROC=completeROC, legends=legends, reduce=reduce, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, xlim = xlim, ylim = ylim, cex.point = cex.point, lwd.curve = lwd.curve, mar = mar, lim.density = lim.density, xlab=xlab, main.density = main.density, col.controlscases = col.controlscases, col.curve = col.curve, col.threshold = col.threshold)
          }
        }
        if(videobar==TRUE){setTxtProgressBar(bar, i)}
        if(!save)  Sys.sleep(tpause)
      })

    }else{
      sapply(1:length(fpr), function(i){
        FPR <- fpr[i]
        plot.buildROC(x, FPR=FPR, h=h, build.process=TRUE, histogram=histogram, breaks=breaks, completeROC=completeROC, legends=legends, reduce=reduce, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, xlim = xlim, ylim = ylim, cex.point = cex.point, lwd.curve = lwd.curve, mar = mar, lim.density = lim.density, xlab=xlab, main.density = main.density, col.controlscases = col.controlscases, col.curve = col.curve, col.threshold = col.threshold)
        if(videobar==TRUE){setTxtProgressBar(bar, i)}
        if(!save)  Sys.sleep(tpause)
      })
    }
    if(videobar==TRUE){close(bar)}

  }

  if(save){
    saveGIF(movie(x, fpr=fpr), movie.name = file, img.name = "Rplot", interval=interval, ani.width = ani.width, ani.height = ani.height, ...)
  }else{
    movie(x, fpr=fpr)
  }

}


movieROC.biroc <- function(x, fpr = NULL, border = TRUE, col = 'blue', completeROC = FALSE, videobar = TRUE, file = "animation1.gif", save = TRUE, legends = FALSE, tpause = 1, interval = 0.2, ani.width = 900, ani.height = 500, cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, alpha.points = 1, alpha.contour = 0.25, lwd.curve = 2, lty.curve = 1, col.curve = 'black', xlab = "X1", ylab = "X2", main.density = "Density functions", lf = NULL, ...){

  col <- col; col.curve <- col.curve; border <- border; completeROC <- completeROC; legends <- legends; videobar <- videobar; tpause <- tpause; interval <- interval; cex <- cex; cex.lab <- cex.lab; cex.axis <- cex.axis; cex.main <- cex.main; alpha.points <- alpha.points; alpha.contour <- alpha.contour; lwd.curve <- lwd.curve; lty.curve <- lty.curve; xlab <- xlab; ylab <- ylab; lf <- lf;  main.density <- main.density

  movie <- function(x, fpr){

    t <- x$t;

    if(is.null(fpr)){
      if(length(t) < 150 ) fpr <- t else fpr <- seq(min(t), max(t), length.out=100)
    }

    B <- length(fpr)

    if(videobar==TRUE){
      cat("\nProgress bar: Construction of GIF with ", B, " thresholds. \n", sep = "")
      bar <- txtProgressBar(min = 0, max = B, style = 3)
    }

    par(mfrow=c(1,2))
    par(mar=c(5.1,5.1,4.1,2.1))
    plot.buildROC(x, col = col, col.curve = "white", border = border, completeROC = completeROC, legends = legends, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main, alpha.points = alpha.points, alpha.contour = alpha.contour, lwd.curve = lwd.curve, lty.curve = lty.curve, xlab = xlab, ylab = ylab, main.density = main.density, lf = lf)

    sapply(1:length(fpr), function(i){
      FPR <- fpr[i]
      plot.buildROC(x, FPR = FPR, build.process = TRUE, col = col, border = border, completeROC = completeROC, legends = legends, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main, alpha.points = alpha.points, alpha.contour = alpha.contour, lwd.curve = lwd.curve, lty.curve = lty.curve, col.curve = col.curve, xlab = xlab, ylab = ylab, main.density = main.density, lf = lf)
      if(videobar==TRUE){setTxtProgressBar(bar, i)}
      if(!save)  Sys.sleep(tpause)
    })
    if(videobar==TRUE){close(bar)}

  }

  if(save){
    saveGIF(movie(x, fpr=fpr), movie.name = file, img.name = "Rplot", interval=interval, ani.width = ani.width, ani.height = ani.height, ...)
  }else{
    movie(x, fpr=fpr)
  }

}


movieROC.multiroc <- function(x, fpr = NULL, display.method = c("PCA", "OV"), displayOV = c(1,2), border = TRUE, col = 'blue', completeROC = FALSE, videobar = TRUE, file = "animation1.gif", save = TRUE, legends = FALSE, tpause = 1, interval = 0.2, ani.width = 900, ani.height = 500, cex = 0.8, cex.lab = 1.5, cex.axis = 1.5, cex.main = 2, alpha.points = 1, alpha.contour = 0.25, lwd.curve = 2, lty.curve = 1, col.curve = 'black', xlab = NULL, ylab = NULL, main.density = "Density functions", lf = NULL, ...){

  display.method <- match.arg(display.method)
  col <- col; col.curve <- col.curve; border <- border; completeROC <- completeROC; legends <- legends; videobar <- videobar; tpause <- tpause; interval <- interval; cex <- cex; cex.lab <- cex.lab; cex.axis <- cex.axis; cex.main <- cex.main; alpha.points <- alpha.points; alpha.contour <- alpha.contour; lwd.curve <- lwd.curve; lty.curve <- lty.curve; xlab <- xlab; ylab <- ylab; lf <- lf;  main.density <- main.density

  movie <- function(x, fpr){

    t <- x$t;

    if(is.null(fpr)){
      if(length(t) < 150 ) fpr <- t else fpr <- seq(min(t), max(t), length.out=100)
    }

    B <- length(fpr)

    if(videobar==TRUE){
      cat("\nProgress bar: Construction of GIF with ", B, " thresholds. \n", sep = "")
      bar <- txtProgressBar(min = 0, max = B, style = 3)
    }

    par(mfrow=c(1,2))
    par(mar=c(5.1,5.1,4.1,2.1))
    plot.buildROC(x, display.method = display.method, displayOV = displayOV, col = col, col.curve = "white", border = border, completeROC = completeROC, legends = legends, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main, alpha.points = alpha.points, alpha.contour = alpha.contour, lwd.curve = lwd.curve, lty.curve = lty.curve, xlab = xlab, ylab = ylab, main.density = main.density)

        sapply(1:length(fpr), function(i){
      FPR <- fpr[i]
      plot.buildROC(x, display.method = display.method, displayOV = displayOV, FPR = FPR, build.process = TRUE, col = col, border = border, completeROC = completeROC, legends = legends, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main, alpha.points = alpha.points, alpha.contour = alpha.contour, lwd.curve = lwd.curve, lty.curve = lty.curve, col.curve = col.curve, xlab = xlab, ylab = ylab, main.density = main.density)
      if(videobar==TRUE){setTxtProgressBar(bar, i)}
      if(!save)  Sys.sleep(tpause)
    })
    if(videobar==TRUE){close(bar)}

  }

  if(save){
    saveGIF(movie(x, fpr=fpr), movie.name = file, img.name = "Rplot", interval=interval, ani.width = ani.width, ani.height = ani.height, ...)
  }else{
    movie(x, fpr=fpr)
  }

}

