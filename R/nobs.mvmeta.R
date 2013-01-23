###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2012-2013
#
`nobs.mvmeta` <-
function (object, ...) {
#
################################################################################
# EXTRACTS THE NUMBER OF OBSERVATIONS USED FOR FITTING. USED BY BIC
#
   return(object$df$nobs)
#
}