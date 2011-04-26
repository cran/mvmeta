fitted.mvmeta <-
function(object, format=c("matrix","list"), na.action, ...) {

	format <- match.arg(format,c("matrix","list"))
	if(missing(na.action)) na.action <- object$na.action
	na.action <- match.arg(na.action,c("na.omit","na.exclude",
		"na.fail","na.pass"))

	# RE-CREATE OBJECTS
	# HERE nalist INCLUDES ALL THE OBS, ALSO MISSING
	# MISSING VALUES IN X ARE EXCLUDED IN A SECOND STAGE
	nalist <- lapply(object$y, function(x) rep(TRUE,length(x)))
	kXlist <- with(object,kXlistmk(X,cen,nalist,length(y),dim$k))

	# COMPUTE FITTED
	fitted <- lapply(kXlist,function(kX) {
		fit <- as.numeric(kX%*%object$beta)
		names(fit) <- object$lab$klab
		return(fit)})

	# FORMAT
	if(format=="matrix") {
		fitted <- rbindlist(fitted)
		rownames(fitted) <- object$lab$mlab
	} else names(fitted) <- object$lab$mlab

	# HANDLE MISSING
	if(!is.null(object$X)) {Xna <- !is.na(object$X)
	} else Xna <- as.matrix(rep(TRUE,length(object$y)))
	nastudy <- rowSums(Xna)!=0
	if(na.action=="na.fail" && !all(Xna)) {
		stop("missing values in 'X'")
	}
	if(na.action=="na.omit") {
		if(format=="matrix") { fitted <- fitted[nastudy,,drop=FALSE]
		} else fitted <- fitted[nastudy]
	}

	return(fitted)
}

