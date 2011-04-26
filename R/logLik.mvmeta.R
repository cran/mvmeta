logLik.mvmeta <-
function(object, ...) {

	val <- object$logLik

	attr(val, "nall") <- object$df$nobs
	attr(val, "nobs") <- object$df$nobs - 
		(object$method=="reml")*object$df$fixed
	attr(val, "df") <- object$df$fixed + object$df$random
	attr(val, "fixed") <- object$df$fixed
	attr(val, "random") <- object$df$random
	class(val) <- "logLik"
	return(val)
}

