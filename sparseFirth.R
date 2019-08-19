###NOTES###
# theta = (beta, c) is the full parameter matrix
# X = [X, D]  is design matrix plus group dummies

spInfoLogit <- function(theta, y, X){
  # Compute the estimated Fisher Information
  # to be used for the BR penalty
  XB <- drop(X %*% theta) #not sparse
  W <- plogis(XB) * plogis(XB, lower.tail=F)
  WX <- sqrt(W)*X
  return(crossprod(WX))
}      
spBrLogit <- function(theta, y, X){
  #Log likelihood function for FFE
  XB <- drop(X %*% theta) #not sparse
  D <- determinant(spInfoLogit(theta, y, X))
  br <- sum(y*(XB) - log(1 + exp(XB))) + .5*D$mod[1]
  return(-br)
}



hat.spam <- function(WX){ #Spam is faster, but doesn't always work
  # Return the diagonals of the logit Hat matrix
  # h= diag(W^{1/2}X(X'WX)^{-1}(W^{1/2}X)).
  # If WX is sqrt(W)X, then this problem is
  # WX ((WX)'WX)^{-1} (WX)' OR
  # WX A^{-1} WX' 
  # WX U^{-1} U'^{-1} WX'  by Cholesky
  # Z'Z
  # Solve U'Z = WX by forward substitution and we're there
  
  WX <- as.spam.dgCMatrix(WX)
  A <- spam::crossprod(WX)
  U <- spam::chol(A)
  Z <- spam::forwardsolve(t(U), t(WX[,ordering(U)]))
  h <- colSums(Z^2)
  return(h)
}

hat.Matrix <- function(WX){ #Matrix always works, but is slower
  A <- Matrix::crossprod(WX)
  U <- Matrix::chol(A, pivot=T)
  Z <- forwardsolve(t(U), t(WX[, attr(U, "pivot")]))
  h <- colSums(Z^2)
  return(h)
}


spGrBrLogit <- function(theta, y, X){
  XB <- drop(X %*% theta) #not sparse
  p1 <- plogis(XB)
  p0 <- plogis(XB, lower.tail=F)
  W <-  p1 * p0
  WX <- sqrt(W)*X
  
  #Try spam first. If that failures, use slower Matrix version
  h <- tryCatch(hat.spam(WX), error=function(x){hat.Matrix(WX)})  
  
  #If the both fail, the whole gradient is a failure
  if(class(h)!="numeric"){stop("error in computing h")}
  
  #Score function 
  Score <- (y- p1)*X*(1+h/2) + (p0-y)*X*(h/2)
  return(-colSums(Score))
}


brLogit.OPG <- function(theta, y, X){
  XB <- drop(X %*% theta) 
  p <- plogis(XB)
  p1 <- plogis(XB, lower.tail=F)
  W <-  p* p1
  WX <- sqrt(W)*X
  h <- tryCatch(hat.spam(WX), error=function(x){hat.Matrix(WX)})  
  if(class(h)!="numeric"){stop("error in computing h")}
  Score <- (y- p)*X*(1+h/2) + (p1-y)*X*(h/2)
  return(-Score) #gradient for each observation for the BHHH estimator
}




sp.ffe <- function(formula, data, subset, na.action=na.omit, start, optim.control=list() ){
  cl <- match.call()
  
  ## make the model frame
  mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
  mf <- cl[c(1L, mf)]
  mf$formula <- formula
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  X <- sparse.model.matrix(formula,data=mf)
  y <- model.response(mf, "any")
  
  if(missing(start)){
    start <- rnorm(ncol(X))*.01
  }
  names(start) <- colnames(X)
  est <- optim(start,spBrLogit,
               gr=spGrBrLogit, 
               # hessian=TRUE,
               X=X, y=y, method="BFGS",
               control=optim.control)
  
  A <- spInfoLogit(est$par, y, X)
  V <- solve(A) #I(theta) = (G'G)^{-1}, solve works well here
  se <- sqrt(diag(V)) 
  t <- est$par/se
  p <- pnorm(abs(t), lower=F)*2
  fitted.values <- plogis(drop(X %*% est$par))
  
  out <- list(coefficients=est$par, 
              se=se, 
              logLik = -est$value,
              vcov=V,
              call=cl,
              X = X, 
              y = y,
              fitted.values = fitted.values,
              convergence=est$convergence,
              coef.table=cbind(est$par, se, t, p))
  
  return(out)
}





spLogit <- function(theta, y, X){
  #Log likelihood function for ordinary logit
  XB <- drop(X %*% theta) #not sparse
  br <- sum(y*(XB) - log(1 + exp(XB))) 
  return(-br)
}


spGrLogit <- function(theta, y, X){
  XB <- drop(X %*% theta) #not sparse
  p1 <- plogis(XB)
  p0 <- plogis(XB, lower.tail=F)
  
  #Score function 
  Score <- (y*p0 -(1-y)*p1)*X
  return(-colSums(Score))
}


brLogit.OPG <- function(theta, y, X){
  XB <- drop(X %*% theta) #not sparse
  p1 <- plogis(XB)
  p0 <- plogis(XB, lower.tail=F)
  
  #Score function 
  Score <- (y*p0 -(1-y)*p1)*X
  return(-Score) #gradient for each observation for the BHHH estimator
}



sp.logit <- function(formula, data, subset, na.action=na.omit, start, optim.control=list() ){
  cl <- match.call()
  
  ## make the model frame
  mf <- match(c("data", "subset", "na.action"), names(cl), 0L)
  mf <- cl[c(1L, mf)]
  mf$formula <- formula
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  X <- sparse.model.matrix(formula,data=mf)
  y <- model.response(mf, "any")
  
  if(missing(start)){
    start <- rnorm(ncol(X))*.01
  }
  names(start) <- colnames(X)
  est <- optim(start,spLogit,
               gr=spGrLogit, 
               # hessian=TRUE,
               X=X, y=y, method="BFGS",
               control=optim.control)
  
  A <- spInfoLogit(est$par, y, X)
  V <- solve(A) #I(theta) = (G'G)^{-1}, solve works well here
  se <- sqrt(diag(V)) 
  t <- est$par/se
  p <- pnorm(abs(t), lower=F)*2
  fitted.values <- plogis(drop(X %*% est$par))
  
  out <- list(coefficients=est$par, 
              se=se, 
              logLik = -est$value,
              vcov=V,
              call=cl,
              X = X, 
              y = y,
              fitted.values = fitted.values,
              convergence=est$convergence,
              coef.table=cbind(est$par, se, t, p))
  
  return(out)
}



