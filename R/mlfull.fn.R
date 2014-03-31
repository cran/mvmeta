###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
mlfull.fn <-
function(par, Xlist, ylist, Slist, nalist, k, m, p, nall, bscov, ctrl) {
#
################################################################################
#
  # COMPUTE Psi FROM PARAMETERS DEPENDING ON STRUCTURE AND PARAMETERIZATION
  Psi <- par2Psi(par[-seq(k*p)],k,bscov,ctrl)
#  
  # EXTRACT coef
  coef <- par[seq(k*p)]
#
  # FIT BY GLS (ONLY USED FOR QUANTITIES, NOT COEF)
  gls <- glsfit(Xlist,ylist,Slist,nalist,Psi,onlycoef=FALSE)
#
  # LIKELIHOOD FUNCTION
  # CONSTANT PART
  pconst <- -0.5*nall*log(2*pi)
  # RESIDUAL COMPONENT
  pres <- -0.5*(crossprod(gls$invtUy-gls$invtUX%*%coef))
  # DETERMINANT COMPONENT
  pdet <- -sum(sapply(gls$Ulist,function(U) sum(log(diag(U)))))
#
  # RETURN
  as.numeric(pconst + pdet + pres)
}
