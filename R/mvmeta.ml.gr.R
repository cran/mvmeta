###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`mvmeta.ml.gr` <-
function(par, Xlist, ylist, Slist, nalist, k, m, nall) {
#
################################################################################
#
  # PARAMETERIZATION OF theta AS THE LOWER TRIANGULAR COMPONENTS OF
  #  THE SQUARE ROOT OF Psi, THROUGH CHOLESKY-DECOMPOSITION
  L <- diag(0,k)
  L[lower.tri(L, diag = TRUE)] <- par
  U <- t(L)
  Psi <- crossprod(U)
# 
  # FIT BY GLS
  gls <- .gls(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#
  # COMPUTE QUANTITIES
  invSigmalist <- lapply(gls$invUlist,tcrossprod)
  reslist <- mapply(function(X,y) as.numeric(y-X%*%gls$coef),
    Xlist,ylist,SIMPLIFY=FALSE)
  ind1 <- rep(1:k,k:1)
  ind2 <- unlist(sapply(1:k,seq,to=k))
#
  # COMPUTE THE MATRIX DERIVATIVES OF EACH PARAMETER
  grad <- .gradchol.ml(par,U,ind1,ind2,invSigmalist,reslist,nalist,k,m)
#
  return(grad)
}

