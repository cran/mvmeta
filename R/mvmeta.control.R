###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012
#
mvmeta.control <- 
  function (fnscale=-1, ...) {
#
################################################################################
# SET CONTROL PARAMETERS FOR MODEL FITTING, WITH SPECIFIC DEFAULT VALUES
	return(list(fnscale=fnscale, ...))
}