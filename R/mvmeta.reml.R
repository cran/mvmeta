mvmeta.reml <-
function(par, ylist, Slist, kXlist, nalist, nobs, k) {

	# PARAMETERIZATION OF theta AS THE LOWER TRIANGULAR COMPONENTS OF
	#	THE SQUARE ROOT OF Psi, THROUGH CHOLESKY-DECOMPOSITION
	Psi <- matrix(0,k,k)
	Psi[lower.tri(Psi,diag=TRUE)] <- par
	Psi <- tcrossprod(Psi)

	# GET THE ESTIMATE OF coef CONDITIONAL ON Psi_theta
	Sigmalist <- mapply(function(S,na) S+Psi[na,na,drop=FALSE],
		Slist,nalist,SIMPLIFY=FALSE)
	Ulist <- lapply(Sigmalist,chol)
	invUlist <- lapply(Ulist,function(U) backsolve(U,diag(ncol(U))))
	invtUXlist <- mapply(function(invU,kX) crossprod(invU,kX),
		invUlist,kXlist,SIMPLIFY=FALSE)
	invtUylist <- mapply(function(invU,y) crossprod(invU,y),
		invUlist,ylist,SIMPLIFY=FALSE)
	coef <- qr.solve(rbindlist(invtUXlist),rbindlist(invtUylist))

	# LIKELIHOOD FUNCTION
	# CONSTANT PART
	const <- -0.5*(nobs-length(coef))*log(2*pi)
	# I GUESS IN STATA:
	#const <- -0.5*(length(ylist)-1)*ncol(Psi)*log(2*pi)
	# RESIDUAL COMPONENT
	res <- -0.5*sum(mapply(function(invtUy,invtUX) {
		crossprod(invtUy-invtUX%*%coef)},invtUylist,invtUXlist))
	# DETERMINANT COMPONENTS
	det1 <- -sum(sapply(Ulist,function(U) sum(log(diag(U)))))
	tXMXtot <- sumlist(lapply(invtUXlist,function(x)crossprod(x)))
	det2 <- -sum(log(diag(chol(tXMXtot))))

	logLik <- const + det1 + det2 + res

	# RETURN MINUS THE VALUE: optim MINIMIZES
	return(logLik)
}

