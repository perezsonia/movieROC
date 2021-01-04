buildROC <- function(X, D, ...) {
  UseMethod("buildROC")
}


buildROC.univarDens <- function(X, D, videobar=FALSE, ...){

  levels.names <- levels(as.factor(D))
  controls <- split(X,D)[[levels.names[1]]]; cases <- split(X,D)[[levels.names[2]]]

  cut.off <- round(c(seq(-4,-2.1,0.09),seq(-2,2,0.06),seq(2.1,5.5,0.09)), 2)
  B <- length(cut.off)

  if(videobar==TRUE){
    cat("\nProgress bar: Construction of GIF with ", B, " thresholds. \n", sep = "")
    bar <- txtProgressBar(min = 0, max = B, style = 3)
  }

  for(i in 1:length(cut.off)){
    if (i < 10) {name = paste('000',i,'plot.png',sep='')}
    if (i < 100 && i >= 10) {name = paste('00',i,'plot.png', sep='')}
    if (i >= 100) {name = paste('0', i,'plot.png', sep='')}

    C <- cut.off[i]
    plotdensityROC(controls, cases, C, build.process=TRUE)
    # dev.off()
    if(videobar==TRUE){setTxtProgressBar(bar, i)}

    close(bar)
  }

  plotdensityROC(controls, cases)
}