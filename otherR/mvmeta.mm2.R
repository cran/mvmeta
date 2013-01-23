###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`mvmeta.mm2` <-
function(Xlist, ylist, Slist, nalist, k, m, p, nall, control) {
#
################################################################################
#
  # FIT FIXED EFFECTS MODEL
  Psi <- diag(0,k)
  gls <- .gls(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#
  tXWXtot <- .sumlist(lapply(gls$invtUXlist,crossprod))
  invtXWXtot <- chol2inv(chol(tXWXtot))
#
  # COMPUTE Phi
  Philist <- mapply(function(invtUX,invU,X) tcrossprod(invU) - invU%*%invtUX %*%
    invtXWXtot %*% t(invU%*%invtUX),gls$invtUXlist,gls$invUlist,Xlist,
    SIMPLIFY=FALSE)
  Phi <- matrix(0,k,k)
  for(i in seq(m)) Phi[!nalist[[i]],!nalist[[i]]] <- 
    Phi[!nalist[[i]],!nalist[[i]]] + Philist[[i]]
  # COMPUTE R
  Rlist <- mapply(function(invU,y,X) tcrossprod(invU) %*% 
    tcrossprod(y-X%*%gls$coef),gls$invUlist,ylist,Xlist,SIMPLIFY=FALSE)
  R <- matrix(0,k,k)
  for(i in seq(m)) R[!nalist[[i]],!nalist[[i]]] <- 
    R[!nalist[[i]],!nalist[[i]]] + Rlist[[i]]
  # COMPUTE mI
  mI <- diag(colSums(matrix(!unlist(nalist),m,k,byrow=T))-p,k)
#
  # SOLVE THE SYSTEM
  Psi1 <- qr.solve(Phi,R-mI)
#
  # CREATE Psi
  dim(Psi1) <- c(k,k)
  Psi <- (Psi1+t(Psi1))/2
#
  # FORCE POSITIVE SEMI-DEFINITENESS
  eig <- eigen(Psi)
  Psi <- eig$vectors %*% diag(pmax(eig$values,0),k) %*% t(eig$vectors)
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
    fitted.values=fitted,df.residual=nall-rank-length(par),rank=rank,logLik=NA,
    control=control)
#
  return(fit)
}
