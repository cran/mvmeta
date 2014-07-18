###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2014
#
mvmeta.control <- 
function(optim=list(), showiter=FALSE, maxiter=100, initPsi=NULL, Psifix=NULL,
  Psicor=0, Scor=0, augment=FALSE, augvar=10^4, igls.iter=10, vc.adj=TRUE,
  reltol=sqrt(.Machine$double.eps), set.negeigen=sqrt(.Machine$double.eps)) {
#
################################################################################
# SET CONTROL PARAMETERS FOR MODEL FITTING, WITH SPECIFIC DEFAULT VALUES
#
  # OPTIM:
  optim <- modifyList(optim,list(fnscale=-1,maxit=maxiter,reltol=reltol))
  if(showiter) {
    optim$trace <- 6
    optim$REPORT <- 1
  }
#
  if(igls.iter<1) stop("'igls.iter' in the control list must be positive")
#
  # RETURN
	list(optim=optim,showiter=showiter,maxiter=maxiter,initPsi=initPsi,
    Psifix=Psifix,Psicor=Psicor,Scor=Scor,augment=augment,augvar=augvar,
    igls.iter=igls.iter,vc.adj=vc.adj,reltol=reltol,set.negeigen=set.negeigen)
}
