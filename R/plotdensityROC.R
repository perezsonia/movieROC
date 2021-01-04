plotdensityROC <- function(X, D, ...) {
  UseMethod("plotdensityROC")
}

plotdensityROC.default <- function(X, D, C=NULL, side=c('right', 'left'), build.process=FALSE, h=1, completeROC=TRUE, legends=FALSE, rel.tol=1e-3, par.specify = FALSE, cex.lab = 1.5, cex.axis = 1.25, cex.main = 1.75, lwd = 2, ...){

  levels.names <- levels(as.factor(D))
  controls <- split(X,D)[[levels.names[1]]]; cases <- split(X,D)[[levels.names[2]]]
  side <- match.arg(side)

  plot.density <- function(controls, cases, h){
    par(mar = c(4.1, 4.6, 3.1, 1.1))
    plot(density(controls, adjust=h)$x, density(controls, adjust=h)$y, xlim=c(min(c(density(controls, adjust=h)$x, density(cases, adjust=h)$x)), max(c(density(controls, adjust=h)$x, density(cases, adjust=h)$x))), ylim=c(-0.05-max(density(cases, adjust=h)$y), 0.05+max(density(controls, adjust=h)$y)), 'l', col='darkolivegreen3', lwd = lwd, xlab="Marker (x)", ylab="f(x)", main="Density functions", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    lines(density(cases, adjust=h)$x, -density(cases, adjust=h)$y, col='red', lwd = lwd)
    abline(h=0, col='gray', lwd = 1)
    if(legends) legend('topright', inset=0.03, c("Controls", "Cases"), col=c('darkolivegreen3', 'red'), lty=1)
  }

  f.controls <- approxfun(density(controls, adjust=h)$x, density(controls, adjust=h)$y, yleft=0, yright=0)
  f.cases <- approxfun(density(cases, adjust=h)$x, density(cases, adjust=h)$y, yleft=0, yright=0)

  plot.ROC.density <- function(controls, cases, side, build.process, p.ROC.C=0, Se.ROC.C=0, h){
    par(mar = c(4.1, 4.6, 3.1, 1.1))
    mm <- min(c(density(controls, adjust=h)$x,density(cases, adjust=h)$x)); MM <-  max(c(density(controls, adjust=h)$x,density(cases, adjust=h)$x))
    points <- seq(mm, MM, length.out=200)
    integrate.density.p.ROC <- sapply(points, function(C) integrate(f.controls, ifelse(side=='right',C,mm), ifelse(side=='right',Inf,C), rel.tol = rel.tol, subdivisions = 1000, stop.on.error = FALSE)$value)
    integrate.density.Se.ROC <- sapply(points, function(C) integrate(f.cases, ifelse(side=='right',C,mm), ifelse(side=='right',Inf,C), rel.tol = rel.tol, subdivisions = 1000, stop.on.error = FALSE)$value)
    if(build.process){
      plot(integrate.density.p.ROC, integrate.density.Se.ROC, 'l', lwd = lwd, col=ifelse(completeROC,'gray','white'), xlab="1-Specificity", ylab="Sensitivity", main="ROC curve", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
      if(side=='right'){
        lines(integrate.density.p.ROC[integrate.density.p.ROC>=p.ROC.C & integrate.density.Se.ROC>=Se.ROC.C], integrate.density.Se.ROC[integrate.density.p.ROC>=p.ROC.C & integrate.density.Se.ROC>=Se.ROC.C], lwd = 2*lwd)
      }else{
        lines(integrate.density.p.ROC[integrate.density.p.ROC<=p.ROC.C & integrate.density.Se.ROC<=Se.ROC.C], integrate.density.Se.ROC[integrate.density.p.ROC<=p.ROC.C & integrate.density.Se.ROC<=Se.ROC.C], lwd = 2*lwd)
      }
    }else{
      plot(integrate.density.p.ROC, integrate.density.Se.ROC, 'l', lwd = 2*lwd, xlab="1-Specificity", ylab="Sensitivity", main="ROC curve", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    }
    abline(0, 1, col='gray')
    if(side=='right'){
      xaxis.C <- integrate.density.p.ROC[integrate.density.p.ROC>=p.ROC.C & integrate.density.Se.ROC>=Se.ROC.C]
      yaxis.C <- points[integrate.density.p.ROC>=p.ROC.C & integrate.density.Se.ROC>=Se.ROC.C]
    }else{
      xaxis.C <- integrate.density.p.ROC[integrate.density.p.ROC<=p.ROC.C & integrate.density.Se.ROC<=Se.ROC.C]
      yaxis.C <- points[integrate.density.p.ROC<=p.ROC.C & integrate.density.Se.ROC<=Se.ROC.C]
    }
    list(xaxis.C=xaxis.C, yaxis.C=yaxis.C)
  }

  if(!par.specify){par(mfrow = c(1, 2))}

  mm <- min(c(density(controls, adjust=h)$x,density(cases, adjust=h)$x))
  MM <- max(c(density(controls, adjust=h)$x,density(cases, adjust=h)$x))
  ord.round <- 2 - round(log10(MM-mm))

  if(!is.null(C) && C < mm){
    print("C is lower than all marker values")
    C <- NULL
  }else if(!is.null(C) && C > MM){
    print("C is larger than all marker values")
    C <- NULL
  }

  if(is.null(C)){
    plot.density(controls, cases, h=h)
    ROC <- plot.ROC.density(controls, cases, side=side, h=h, build.process=FALSE)
  }else{
    m <- matrix(c(1,3,3,2,4,4), 3, 2)
    layout(m)
    par(mar = c(4, 4, 2, 1))

    plot.density(controls, cases, h=h)
    if(side=='right'){
      polygon(c(C, density(controls, adjust=h)$x[density(controls, adjust=h)$x > C], tail(density(controls, adjust=h)$x[density(controls, adjust=h)$x > C],1)), c(0, density(controls, adjust=h)$y[density(controls, adjust=h)$x > C], tail(density(controls, adjust=h)$y[density(controls, adjust=h)$x > C],1)), col='darkolivegreen3', border=NA)
      polygon(c(C, density(cases, adjust=h)$x[density(cases, adjust=h)$x > C], tail(density(cases, adjust=h)$x[density(cases, adjust=h)$x > C],1)), c(0, -density(cases, adjust=h)$y[density(cases, adjust=h)$x > C], tail(-density(cases, adjust=h)$y[density(cases, adjust=h)$x > C],1)), col='red', border=NA)
    }else{
      polygon(c(C, density(controls, adjust=h)$x[density(controls, adjust=h)$x < C], tail(density(controls, adjust=h)$x[density(controls, adjust=h)$x < C],1)), c(0, density(controls, adjust=h)$y[density(controls, adjust=h)$x < C], tail(density(controls, adjust=h)$y[density(controls, adjust=h)$x < C],1)), col='darkolivegreen3', border=NA)
      polygon(c(C, density(cases, adjust=h)$x[density(cases, adjust=h)$x < C], tail(density(cases, adjust=h)$x[density(cases, adjust=h)$x < C],1)), c(0, -density(cases, adjust=h)$y[density(cases, adjust=h)$x < C], tail(-density(cases, adjust=h)$y[density(cases, adjust=h)$x < C],1)), col='red', border=NA)
    }

    if(legends) legend('bottomright', inset=0.03, c("1 - Sp", "Se"), fill=c('darkolivegreen3', 'red'))
    abline(v=C, col='blue', lty=4)
    text(C, -0.05-max(density(cases, adjust=h)$y), adj=1, round(C, ord.round), cex=1, col='blue')

    p.ROC.C <- integrate(f.controls, ifelse(side=='right',C,mm), ifelse(side=='right',Inf,C), rel.tol = rel.tol, subdivisions = 1000, stop.on.error = FALSE)$value
    Se.ROC.C <- integrate(f.cases, ifelse(side=='right',C,mm), ifelse(side=='right',Inf,C), rel.tol = rel.tol, subdivisions = 1000, stop.on.error = FALSE)$value
    ROC.until.C <- plot.ROC.density(controls, cases, side=side, build.process, p.ROC.C=p.ROC.C, Se.ROC.C=Se.ROC.C, h=h)
    points(p.ROC.C, Se.ROC.C, cex=1.25, pch = 19, col='blue')
    text(p.ROC.C, Se.ROC.C-0.025, round(C, ord.round), adj=0, cex=1, col='blue')

    plot(rep(-0.75,length(controls)), controls, 'p', pch=1, col='green', xlim=c(-2.5,2.5), ylim=c(mm,MM), xaxt='n', xlab="", ylab="Marker values", xaxs='i', cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    points(rep(0.75,length(cases)), cases, col='red')
    boxplot(controls, at=-1.5, col='green', add=TRUE, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    boxplot(cases, at=1.5, col='red', add=TRUE, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    abline(h=C, col='blue', lty=4)
    polygon(c(2.3,2.3,2.5,2.5), c(ifelse(side=='right',C, mm),ifelse(side=='right',MM,C),ifelse(side=='right',MM,C),ifelse(side=='right',C,mm)), col='gray', border=NA)
    if(legends) legend('topright', inset=0.03, c("Controls", "Cases"), col=c('green', 'red'), pch=c(1,1))

    plot(0.5, controls[1], lwd = 2*lwd, col='white', xlim=c(0,1), ylim=c(mm,MM), xlab="1-Specificity", ylab="Marker intervals", main="Classification subsets", cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main)
    polygon(c(ROC.until.C$xaxis.C, rev(ROC.until.C$xaxis.C)), c(ROC.until.C$yaxis.C, rep(ifelse(side=='right',MM,mm),length(ROC.until.C$yaxis.C))), col='gray', border=NA)
    if(legends) legend('topright', inset=0.03, c("Classified as Controls","Classified as Cases"), fill=c('white', 'gray'))
  }

}
