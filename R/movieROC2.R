movieROC2 <- function(X, D, ...) {
  UseMethod("movieROC2")
}

movieROC2.univarDens <- function(X, D, side=c('right','left'), h=1, cut.off=NULL, completeROC=FALSE, legends=FALSE, videobar=FALSE, file="animation1.gif", clean=FALSE, interval=0.2, ani.width = 500, ani.height = 750, ...){

  side <- match.arg(side)

  movie <- function(X, D, cut.off, completeROC, videobar, legends){

    if(is.null(cut.off)){
      range <- max(X) - min(X)
      if(length(unique(X)) < 150 ){
        cut.off <- c(min(X) - range/20, sort(unique(X)), max(X) + range/20)
      }else{
        cut.off <- c(min(X) - range/20, seq(min(X), max(X), length.out=100), max(X) + range/20)
      }
    }

    B <- length(cut.off)

    if(videobar==TRUE){
      cat("\nProgress bar: Construction of GIF with ", B, " thresholds. \n", sep = "")
      bar <- txtProgressBar(min = 0, max = B, style = 3)
    }

    sapply(1:length(cut.off), function(i){
      C <- cut.off[i]
      plotdensityROC(X, D, side=side, C, h=h, build.process=TRUE, completeROC=completeROC, legends=legends)
      if(videobar==TRUE){setTxtProgressBar(bar, i)}
    })
    if(videobar==TRUE){close(bar)}
  }

  saveGIF(movie(X, D, cut.off=cut.off, completeROC=completeROC, legends=legends, videobar=videobar), movie.name = file, img.name = "Rplot", clean=clean, interval=interval, ani.width = ani.width, ani.height = ani.height, ...)

}
