###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`qtest.mvmeta` <-
function(object, ...) {
#
################################################################################
# RUN THE RELATED FIXED-EFFECTS MODEL
# 
  mf <- model.frame(object)
  int <- attr(object$terms,"intercept")==1L
  y <- as.matrix(model.response(mf,"numeric"))
  X <- model.matrix(object) 
  S <- object$S
  nay <- is.na(y)
#
  # TRANSFORM X, y, S AND nay IN LISTS
  dim <- object$dim
  Xlist <- lapply(seq(dim$m),function(i) diag(1,dim$k)[!nay[i,],,drop=FALSE]%x%
    X[i,,drop=FALSE])
  ylist <- lapply(seq(dim$m),function(i) y[i,][!nay[i,]])
  Slist <- lapply(seq(dim$m),function(i) 
    xpndMat(S[i,])[!nay[i,],!nay[i,],drop=FALSE])
  nalist <- lapply(seq(dim$m),function(i) nay[i,])
#
  # GLS 
  Psi <- diag(0,dim$k)
  gls <- .gls(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#
################################################################################
# COMPUTE THE STATS
#
  Q <- c(drop(crossprod(gls$invtUy-gls$invtUX%*%gls$coef)),
    colSums(do.call("rbind",mapply(function(y,S,X,na) {
      comp <- rep(0,dim$k)
      comp[!na] <- as.vector((y-X%*%gls$coef)^2 / diag(S))
      return(comp)},ylist,Slist,Xlist,nalist,SIMPLIFY=FALSE))))
#
  df <- c(with(object$df,nall-fixed),colSums(!nay,na.rm=TRUE)-dim$p)
  pvalue <- sapply(seq(Q),function(i) 1-pchisq(Q[i],df[i]))
  names(Q) <- names(df) <- names(pvalue) <- c(".all",object$lab$k)
#
  qstat <- list(Q=Q,df=df,pvalue=pvalue,residual=object$dim$p-int>0L,k=dim$k)
  class(qstat) <- "qtest.mvmeta"
#
  return(qstat)
}
