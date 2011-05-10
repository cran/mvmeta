qtest.mvmeta <-
function(object, ...) {

	# RE-COMPUTE QUANTITIES
	model <- object$model
	y <- as.matrix(model.response(model,"numeric"))
	y <- lapply(seq(nrow(as.matrix(y))), function(i) as.matrix(y)[i,])
	terms <- terms(model)
	X <- model.matrix(terms,model,object$contrast)
	S <- object$S
	nastudy <- sapply(y,function(x) !all(is.na(x)))
	nastudy[sapply(S,function(x) all(is.na(rowSums(x)+colSums(x))))] <- FALSE
	Xna <- rowSums(!is.na(X))==object$dim$p
	nastudy[!Xna] <- FALSE
	nalist <- mapply(function(y,S) {
		na <- !is.na(y)
		na[unique(which(is.na(S),arr.ind=TRUE))] <- FALSE
		return(na)},y[nastudy],S[nastudy],SIMPLIFY=FALSE)
	ylist <- mapply(function(y,na) y[na],y[nastudy],nalist,SIMPLIFY=FALSE)
	Slist <- mapply(function(S,na) S[na,na,drop=FALSE],
		S[nastudy],nalist,SIMPLIFY=FALSE)
	kXlist <- mapply(function(i,na) {diag(1,object$dim$k)[na,,drop=FALSE]%x%
		X[nastudy,,drop=FALSE][i,,drop=FALSE]},
		seq(length(ylist)),nalist,SIMPLIFY=FALSE)

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

