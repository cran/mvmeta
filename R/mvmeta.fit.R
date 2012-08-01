###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012
#
mvmeta.fit <- 
  function(X, y, S, offset=NULL, method="reml", control=list()) {
#
################################################################################
# PREPARE THE DATA
#
  # SET DIMENSIONS AND MISSING
  y <- as.matrix(y)
  nay <- is.na(y)
  k <- ncol(y)
  m <- nrow(y)
  p <- ncol(X)
  nall <- sum(!nay)
  # STORE NAMES
  nk <- colnames(y)
  if(k>1L && is.null(nk)) nk <- paste("y",seq(k),sep="")
  nm <- rownames(y)
  np <- colnames(X)
#
  # EXCLUDE OFFSET (RECYCLING RULE)
  if(!is.null(offset)) y <- y - offset
#
  # TRANSFORM X, y, S AND nay IN LISTS
  Xlist <- lapply(seq(m),function(i) diag(1,k)[!nay[i,],,drop=FALSE]%x%
    X[i,,drop=FALSE])
  ylist <- lapply(seq(m),function(i) y[i,][!nay[i,]])
  Slist <- lapply(seq(m),function(i) {
    Si <- xpndMat(S[i,])[!nay[i,],!nay[i,],drop=FALSE]
    if(any(is.na(Si))) stop("missing pattern in 'y' and S' is not consistent")
    return(Si)
  })
  nalist <- lapply(seq(m),function(i) nay[i,])
#    
################################################################################
# FIT THE MODEL AND SET OBJECTS
#
  # SELECT THE ESTIMATION METHOD
  control <- do.call("mvmeta.control",control)
  fun <- paste("mvmeta",match.arg(method,c("fixed","ml","reml")),sep=".")
  fit <- do.call(fun,list(Xlist=Xlist,ylist=ylist,Slist=Slist,nalist=nalist,
    k=k,m=m,nall=nall,control=control))
#  
  # DEFINE DIMENSIONS, DF AND LABELS
  fit$dim <- list(k=k,m=m,p=p)
  fit$df <- list(nall=nall,nobs=nall-(method=="reml")*fit$rank,
    df=nall-fit$df.residual,fixed=fit$rank,random=ifelse(method=="fixed",0,
    k*(k+1)/2))
  fit$lab <- list(k=nk,p=np)
#
  # RE-INSERT OFFSET IN y AND FITTED VALUES (RECYCLING RULE)
  # RE-INSTATE DIMENSIONS AND PUT LABELS
  temp <- rep(NA,m*k)
  temp[!t(nay)] <- unlist(fit$fitted.values)
  fit$fitted.values <- matrix(temp,m,k,byrow=TRUE)
  if(!is.null(offset)) {
    y <- y + offset
    fit$fitted.values <- fit$fitted.values+offset
  }
  if(method!="fixed") dimnames(fit$Psi) <- list(nk,nk)
  if(k==1L) {
    names(fit$coefficients) <- np
    dimnames(fit$vcov) <- list(np,np)
    fit$fitted.values <- drop(fit$fitted.values)
    fit$residuals <- drop(y-fit$fitted.values)
    names(fit$residuals) <- names(fit$fitted.values) <- nm
  } else {
    fit$coefficients <- matrix(fit$coefficients,p,k,dimnames=list(np,nk))
    rownames(fit$vcov) <- colnames(fit$vcov) <- 
      paste(rep(nk,each=p),rep(np,k),sep=".")
    fit$residuals <- y-fit$fitted.values
    dimnames(fit$residuals) <- dimnames(fit$fitted.values) <- list(nm,nk)
  }
#
  return(fit)
}
