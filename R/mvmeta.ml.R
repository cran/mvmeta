###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`mvmeta.ml` <-
function(Xlist, ylist, Slist, nalist, k, m, p, nall, control) {
#
################################################################################
#
  # PRODUCE INITIAL VALUES THROUGH IGLS
  Psi <- diag(0.001,k)
  for(i in seq(control$igls.niter)) {
    Psi <- mvmeta.igls(Psi,Xlist,ylist,Slist,nalist,k,m)
  }
#    
  # PARAMETERIZATION OF theta AS THE LOWER TRIANGULAR COMPONENTS OF
  #  THE SQUARE ROOT OF Psi, THROUGH CHOLESKY-DECOMPOSITION
  par <- vechMat(t(chol(Psi)))
#  
  # MAXIMIZE
  opt <- optim(par,mvmeta.ml.fn,mvmeta.ml.gr,Xlist=Xlist,ylist=ylist,
    Slist=Slist,nalist=nalist,k=k,m=m,nall=nall,method="BFGS",
    control=control$optim)
  if(!(converged <- opt$convergence==0L)) {
    warning("convergence not reached after maximum number of iterations")
  }
#
  # Psi: ESTIMATED BETWEEN-STUDY (CO)VARIANCE MATRIX
  Psi <- matrix(0,k,k)
  Psi[lower.tri(Psi,diag=TRUE)] <- opt$par
  Psi <- tcrossprod(Psi)
# 
  # FIT BY GLS
  gls <- .gls(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#
  # COMPUTE (CO)VARIANCE MATRIX OF coef
  qrinvtUX <- qr(gls$invtUX)
  R <- qr.R(qrinvtUX)
  Qty <- qr.qty(qrinvtUX,gls$invtUy)
  vcov <- tcrossprod(backsolve(R,diag(1,ncol(gls$invtUX))))
#
  # COMPUTE RESIDUALS (LATER), FITTED AND RANK
  res <- NULL
  fitted <- lapply(Xlist,"%*%",gls$coef)
  rank <- qrinvtUX$rank
#
  fit <- list(coefficients=gls$coef,vcov=vcov,Psi=Psi,residuals=res,
    fitted.values=fitted,df.residual=nall-rank-length(par),rank=rank,
    logLik=opt$value,converged=converged,control=control)
#
  return(fit)
}

