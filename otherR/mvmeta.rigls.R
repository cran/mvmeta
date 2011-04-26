mvmeta.rigls <-
function(Psi, ylist, Slist, kXlist, nalist, k, m, p) {

#SOMETHING WRONG....

	# COMPUTE beta BY GLS
	Sigmalist <- mapply(function(S,na) S+Psi[na,na,drop=FALSE],
		Slist,nalist,SIMPLIFY=FALSE)
	Ulist <- lapply(Sigmalist,chol)
	invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
	invtUXlist <- mapply(function(invU,kX) crossprod(invU,kX),
		invUlist,kXlist,SIMPLIFY=FALSE)
	invtUylist <- mapply(function(invU,y) crossprod(invU,y),	
		invUlist,ylist,SIMPLIFY=FALSE)
	invtUX <- rbindlist(invtUXlist)
	invtUy <- rbindlist(invtUylist)
	beta <- as.numeric(qr.solve(invtUX,invtUy))

	# CREATE A MATRIX WITH INDICATOR OF (CO)VAR COMPONENTS
	indMat <- xpndMat(seq(k*(k+1)/2))

	# CREATE THE TRANSFORMED RESPONSE AND DESIGN MATRIX, AND THE 
	#	EXPANDED (CO)VARIANCE MATRIX
	Vlist <- mapply(function(na,S) {
		V <- Psi[na,na,drop=FALSE]+S
		return(V%x%V)},nalist,Slist,SIMPLIFY=FALSE)
	flist <- mapply(function(y,S,kX,invU) {
		# ADD THE BIAS CORRECTION TERM
		V1 <- crossprod(crossprod(invU,kX))
		V2 <- chol(V1)
		V3 <- backsolve(V2,diag(ncol(V2)))
		V4 <- tcrossprod(crossprod(kX,V3))
		return(as.numeric(tcrossprod(y-kX%*%beta)+V4)-as.numeric(S))},
		ylist,Slist,kXlist,invUlist,SIMPLIFY=FALSE)
	Zlist <- lapply(nalist,function(na) {
		z <- as.numeric(indMat[na,na,drop=FALSE])
		return(sapply(unique(z),function(x) as.numeric(z==x)))})

	# CREATE TRANFORMED OBJECTS THROUGH CHOLESKY
	Ulist <- lapply(Vlist,chol)
	invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
	invtUZlist <- mapply(function(invU,Z) crossprod(invU,Z),
		invUlist,Zlist,SIMPLIFY=FALSE)
	invtUflist <- mapply(function(invU,f) crossprod(invU,f),
		invUlist,flist,SIMPLIFY=FALSE)
	invtUZ <- rbindlist(invtUZlist)
	invtUf <- rbindlist(invtUflist)

	# ESTIMATE THE COMPONENTS
	theta <- as.numeric(qr.solve(invtUZ,invtUf))
	Psi <- xpndMat(theta)
	# FORCING POSITIVE-DEFINITENESS
	Psi <- expMat(logMat(Psi,positive=TRUE))

	return(Psi)
}

