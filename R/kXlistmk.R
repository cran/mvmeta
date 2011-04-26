kXlistmk <-
function(X,cen,nalist,m,k) {

	# IF X DOES NOT EXIST, 
	if(is.null(X)) {
		kXlist <- mapply(function(i,na) diag(1,k)[na,,drop=FALSE],
			seq(m),nalist,SIMPLIFY=FALSE)
	# IF INSTEAD X IS A MATRIX
	} else {
		# CENTERING
		kXlist <- cbind(1,sweep(as.matrix(X),2,cen))
		# EXPAND
		kXlist <- mapply(function(i,na) 
			diag(1,k)[na,,drop=FALSE]%x%kXlist[i,,drop=FALSE],
			seq(m),nalist,SIMPLIFY=FALSE)
	}
	return(kXlist)
}

