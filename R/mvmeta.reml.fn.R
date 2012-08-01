mvmeta.reml.fn <-
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
  # RESTRICTED LIKELIHOOD FUNCTION
  # CONSTANT PART
  pconst <- -0.5*(nall-length(gls$coef))*log(2*pi)
  # I GUESS IN STATA:
  #const <- -0.5*length(ylist)*ncol(Psi)*log(2*pi)
  # RESIDUAL COMPONENT
  pres <- -0.5*(crossprod(gls$invtUy-gls$invtUX%*%gls$coef))
  # DETERMINANT COMPONENTS
  pdet1 <- -sum(sapply(gls$Ulist,function(U) sum(log(diag(U)))))
  tXMXtot <- .sumlist(lapply(gls$invtUXlist,function(x)crossprod(x)))
  pdet2 <- -sum(log(diag(chol(tXMXtot))))
#
  logLik <- as.numeric(pconst + pdet1 + pdet2 + pres)
#  
  return(logLik)
}
