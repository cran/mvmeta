###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`mvmeta.control` <- 
function (optim=list(), showiter=FALSE, maxiter=100, igls.iter=10,
  reltol=sqrt(.Machine$double.eps),vc.adj=TRUE,
  set.negeigen=sqrt(.Machine$double.eps)) {
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
	return(list(optim=optim,showiter=showiter,maxiter=maxiter,igls.iter=igls.iter,
    reltol=reltol,vc.adj=vc.adj,set.negeigen=set.negeigen))
}
