###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
remlfull.fn <-
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
  # RESTRICTED LIKELIHOOD FUNCTION
  # CONSTANT PART
  pconst <- -0.5*(nall-length(gls$coef))*log(2*pi)
  # RESIDUAL COMPONENT
  pres <- -0.5*(crossprod(gls$invtUy-gls$invtUX%*%coef))
  # DETERMINANT COMPONENTS
  pdet1 <- -sum(sapply(gls$Ulist,function(U) sum(log(diag(U)))))
  tXWXtot <- sumlist(lapply(gls$invtUXlist,crossprod))
  pdet2 <- -sum(log(diag(chol(tXWXtot))))
#
  # RETURN
  as.numeric(pconst + pdet1 + pdet2 + pres)
}
