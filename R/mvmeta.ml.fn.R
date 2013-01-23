###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`mvmeta.ml.fn` <-
function(par, Xlist, ylist, Slist, nalist, k, m, nall) {
#
################################################################################
#
  # PARAMETERIZATION OF theta AS THE LOWER TRIANGULAR COMPONENTS OF
  #  THE SQUARE ROOT OF Psi, THROUGH CHOLESKY-DECOMPOSITION
  Psi <- diag(0,k)
  Psi[lower.tri(Psi,diag=TRUE)] <- par
  Psi <- tcrossprod(Psi)
#  
  # FIT BY GLS
  gls <- .gls(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#  
  # LIKELIHOOD FUNCTION
  # CONSTANT PART
  pconst <- -0.5*nall*log(2*pi)
  # I GUESS IN STATA:
  #const <- -0.5*length(ylist)*ncol(Psi)*log(2*pi)
  # RESIDUAL COMPONENT
  pres <- -0.5*(crossprod(gls$invtUy-gls$invtUX%*%gls$coef))
  # DETERMINANT COMPONENT
  pdet <- -sum(sapply(gls$Ulist,function(U) sum(log(diag(U)))))
#  
  logLik <- as.numeric(pconst + pdet + pres)
#  
  return(logLik)
}
