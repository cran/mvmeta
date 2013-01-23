###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`mvmeta.control` <- 
function (optim=list(), igls.niter=10) {
#
################################################################################
# SET CONTROL PARAMETERS FOR MODEL FITTING, WITH SPECIFIC DEFAULT VALUES
#
  # OPTIM: FORCE MAXIMIZATION
  optim <- modifyList(optim,list(fnscale=-1))
#
	return(list(optim=optim,igls.niter=igls.niter))
}
