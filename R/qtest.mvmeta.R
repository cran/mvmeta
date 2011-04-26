qtest.mvmeta <-
function(object, ...) {

	# RE-COMPUTE QUANTITIES
	nastudy <- sapply(object$y,function(x) !all(is.na(x)))
	nastudy[sapply(object$S,function(x) all(is.na(rowSums(x)+colSums(x))))] <- FALSE
	if(!is.null(object$X)) {
		Xna <- rowSums(!is.na(object$X))==ncol(object$X)
		nastudy[!Xna] <- FALSE
	}

	nalist <- mapply(function(y,S) {
		na <- !is.na(y)
		na[unique(which(is.na(S),arr.ind=TRUE))] <- FALSE
		return(na)},object$y[nastudy],object$S[nastudy],SIMPLIFY=FALSE)
	ylist <- mapply(function(y,na) y[na],object$y[nastudy],nalist,SIMPLIFY=FALSE)
	Slist <- mapply(function(S,na) S[na,na,drop=FALSE],
		object$S[nastudy],nalist,SIMPLIFY=FALSE)
	kXlist <- kXlistmk(object$X[nastudy,,drop=FALSE],object$cen,nalist,
		length(ylist),object$dim$k)

	# COMPUTE beta WITH FIXED EFFECTS MODEL
	# OBTAIN beta BY GLS
	Ulist <- lapply(Slist,chol)
	invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
	invtUXlist <- mapply(function(invU,kX) crossprod(invU,kX),
		invUlist,kXlist,SIMPLIFY=FALSE)
	invtUylist <- mapply(function(invU,y) crossprod(invU,y),
		invUlist,ylist,SIMPLIFY=FALSE)
	invtUX <- rbindlist(invtUXlist)
	invtUy <- rbindlist(invtUylist)
	beta <- as.numeric(qr.solve(invtUX,invtUy))

	# GENERATE STAT, DF AND P-VALUE
	Q <- as.numeric(crossprod(invtUy-invtUX%*%beta))
	df <- object$df$nobs-object$df$fixed
	pvalue <- 1-pchisq(Q,df)
	residual <- !is.null(object$X)
	k <- object$dim$k

	res <- list(Q=Q,df=df,pvalue=pvalue,residual=residual,k=k)
	class(res) <- "qtest.mvmeta"
	return(res)
}

