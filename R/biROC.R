biROC <- function(X, D, ...) {
  UseMethod("biROC")
}

biROC.default <- function(X, D, method = c("lrm", "fixedLinear", "fixedQuadratic", "dynamicEmpirical", "dynamicMeisner"), formula = "D~poly(X.1,2)+poly(X.2,2)+I(X.1*X.2)", stepModel = TRUE, coefLinear = c(1,1), coefQuadratic = c(1,1,0,1,1), K = 201, alpha = 0.5, approxh = 0.5, multiplier = 2, ...){

  if(ncol(X)!=2){
    stop("X should be a matrix with two columns.")
  }
  X <- na.omit(X)

  levels.names <- levels(as.factor(D))
  if(length(levels.names) != 2){
    stop("D should be a vector with two different levels.")
  }
  controls <- matrix(split(X,D)[[levels.names[1]]], nrow = sum(as.factor(D)==levels.names[1]))
  cases <- matrix(split(X,D)[[levels.names[2]]], nrow = sum(as.factor(D)==levels.names[2]))
  D <- ifelse(as.factor(D)==levels.names[1], 0, 1); levels <- levels(as.factor(D))

  for(i in 1:2){
    assign(paste("X.",i, sep=""), X[,i])
  }

  method <- match.arg(method)

  if(method == "lrm"){

    if(is.null(formula)) stop("formula should be correct formula to run glm as a character.")
    model <- glm(as.formula(formula), family=binomial(link='logit'))
    if(stepModel) invisible(capture.output(model <- step(model)))
    Z <- suppressWarnings(predict(model))
    roc <- gROC(Z, D, side='right')

    results <- list(levels = levels.names, controls = controls, cases = cases, X = X, D = D, method = method, formula = formula, model = model, coefModel = model$coefficients, stepModel = stepModel, Z = Z, t = roc$t, roc = roc$roc, c = roc$c, auc = roc$auc)

  }else if(method == "fixedLinear"){

    if(length(coefLinear)!=2) stop("coefLinear should be a vector with two elements indicating the coefficients of each marker in a linear combination coef[1]*X.1 + coef[2]*X.2.")
    Z <- c(tcrossprod(X, matrix(coefLinear, nrow = 1)))
    roc <- gROC(Z, D, side='right')

    results <- list(levels = levels.names, controls = controls, cases = cases, X = X, D = D, method = method, coefLinear = coefLinear, Z = Z, t = roc$t, roc = roc$roc, c = roc$c, auc = roc$auc)

  }else if(method == "fixedQuadratic"){

    if(length(coefQuadratic)!=5) stop("coefQuadratic should be a vector with five elements indicating the coefficients of each marker in a quadratic combination coef[1]*X.1 + coef[2]*X.2 + coef[3]*(X.1*X.2) + coef[4]*(X.1)^2 + coef[5]*(X.2)^2.")

    x <- X[,1]; y <- X[,2]
    Z <- c(tcrossprod(cbind(x,y,x*y,x^2,y^2), matrix(coefQuadratic, nrow = 1)))
    roc <- gROC(Z, D, side='right')

    results <- list(levels = levels.names, controls = controls, cases = cases, X = X, D = D, method = method, coefQuadratic = coefQuadratic, Z = Z, t = roc$t, roc = roc$roc, c = roc$c, auc = roc$auc)

  }else if(method %in% c("dynamicMeisner", "dynamicEmpirical")){

    m <- nrow(controls)

    t <- seq(0,1,1/m)
    N <- nrow(X)
    t <- t[-1]

    x <- X[,1]; y <- X[,2]
    lx <- (max(x)-min(x))/20; ly <- (max(y)-min(y))/20
    xx <- seq(min(x)-lx,max(x)+lx, length.out=100)
    yy <- seq(min(y)-ly,max(y)+ly, length.out=100)

    if(method == "dynamicMeisner"){

      cat("Progress bar: Estimation of the optimal linear combination for FPR using maxTPR package\n"); flush.console()
      bar <- txtProgressBar(min = 0, max = t[m], style = 3)
      CoefTable <- lapply(t, function(ti){
        setTxtProgressBar(bar, ti)
        invisible(capture.output(coef <- maxTPR::maxTPR(as.data.frame(cbind(D, X)), tval = ti, alpha = alpha, approxh = approxh, multiplier = multiplier)$sTPRrslt[c("coef1","coef2")]))
        biroc <- biROC(X, D, method = "fixedLinear", coefLinear = coef)
        indext <- which(abs(biroc$t - ti) < .Machine$double.eps)
        c <- biroc$c[indext]
        Se <- biroc$roc[indext]
        Z <- biroc$Z
        lf <- max(min(sort(unique(Z))[-1]-sort(unique(Z))[-length(unique(Z))]), .Machine$double.eps)
        f <- outer(xx,yy,function(x,y) c(tcrossprod(cbind(x,y), matrix(biroc$coefLinear, nrow = 1))))
        list(coef = coef, c = c, t = ti, roc = Se, Z = Z, f = f, lf = lf)
      })
      close(bar)

    }else{  # method == "dynamicEmpirical"

      alpha <- seq(-1, 1, length = K)
      rate <- 1/alpha
      G <- cbind(c(rep(1,2*K), rep(-1,2*K)), rep(c(alpha, rate),2))

      cat("Progress bar: Estimation of the optimal linear combination for FPR using the empirical estimator\n"); flush.console()
      bar <- txtProgressBar(min = 0, max = t[m], style = 3)
      CoefTable <- lapply(t, function(ti){
        setTxtProgressBar(bar, ti)
        C <- sapply(1:(4*K), function(j){
          Z <- G[j,1] * X[,1] + G[j,2] * X[,2]
          groc <- gROC(X = Z, D = D)
          index.opt <- which.min(groc$c[groc$t <= ti])
          c(tpr = groc$roc[index.opt], c = groc$c[index.opt], g = G[j,])
        })
        index.maxTPR <- which.max(C[1,])
        coef <- C[3:4,index.maxTPR]
        biroc <- biROC(X, D, method = "fixedLinear", coefLinear = coef)
        Z <- biroc$Z
        lf <- max(min(sort(unique(Z))[-1]-sort(unique(Z))[-length(unique(Z))]), .Machine$double.eps)
        f <- outer(xx,yy,function(x,y) c(tcrossprod(cbind(x,y), matrix(biroc$coefLinear, nrow = 1))))
        list(coef = coef, c = C[2,index.maxTPR], t = ti, roc = C[1,index.maxTPR], Z = Z, f = f, lf = lf)
      })
      close(bar)

    }

    coef <- sapply(1:m, function(i){CoefTable[[i]]$coef})
    c <- sapply(1:m, function(i){CoefTable[[i]]$c})
    roc <- sapply(1:m, function(i){CoefTable[[i]]$roc})
    Z <- sapply(1:m, function(i){CoefTable[[i]]$Z})
    auc <- as.numeric(sum(roc[-m]*(t[-1] - t[-m])) + t[1]*roc[1])

    results <- list(levels = levels.names, controls = controls, cases = cases, X = X, D = D, method = method, coefLinear = coef, Z = Z, CoefTable = CoefTable, t = t, roc = roc, c = c, auc = auc)

  }

  attr(results, 'class') <- 'biroc'

  return(results)

}

