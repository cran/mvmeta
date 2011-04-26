mvmeta.igls <-
function(Psi, ylist, Slist, kXlist, nalist, k, m, p) {

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
	npar <- k*(k+1)/2
	indMat <- xpndMat(seq(npar))

	# CREATE THE TRANSFORMED OBJECTS IN IGLS STRUCTURE
	# EXPANDED (CO)VARIANCE MATRIX
	Sigmalist <- mapply(function(na,S) {
		V <- Psi[na,na,drop=FALSE]+S
		return(V%x%V)},nalist,Slist,SIMPLIFY=FALSE)
	# RESPONSE VECTORS WITH RESIDUALS MINUS THE WITHIN (CO)VARIANCE
	#	COMPONENTS, CONSIDERED FIXED
	flist <- mapply(function(y,S,kX) {
		return(as.numeric(tcrossprod(y-kX%*%beta))-as.numeric(S))},
		ylist,Slist,kXlist,SIMPLIFY=FALSE)
	# DESIGN MATRIX MAPPING THE PARAMETERS TO BE ESTIMATED
	#	IT AUTOMATICALLY CREATES 0 COLUMNS FOR MISSING OBSERVATIONS
	Zlist <- lapply(nalist,function(na) {
		z <- as.numeric(indMat[na,na,drop=FALSE])
		Z <- lapply(seq(npar),function(x) as.numeric(z==x))
		return(cbindlist(Z))})

	# CREATE TRANFORMED OBJECTS FOR WEIGHTED LEAST-SQUARE THROUGH CHOLESKY
	Ulist <- lapply(Sigmalist,chol)
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

