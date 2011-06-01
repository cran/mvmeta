mvmeta.ml <-
function(par, ylist, Slist, kXlist, nalist, nobs, k) {

	# PARAMETERIZATION OF theta AS THE LOWER TRIANGULAR COMPONENTS OF
	#	THE SQUARE ROOT OF Psi, THROUGH CHOLESKY-DECOMPOSITION
	Psi <- matrix(0,k,k)
	Psi[lower.tri(Psi,diag=TRUE)] <- par
	Psi <- tcrossprod(Psi)

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

	# LIKELIHOOD FUNCTION
	# CONSTANT PART
	const <- -0.5*nobs*log(2*pi)
	# I GUESS IN STATA:
	#const <- -0.5*length(ylist)*ncol(Psi)*log(2*pi)
	# RESIDUAL COMPONENT
	res <- -0.5*crossprod(invtUy-invtUX%*%beta)
	# DETERMINANT COMPONENT
	det <- -sum(sapply(Ulist,function(U) sum(log(diag(U)))))

	logLik <- const + det + res

	# RETURN MINUS THE VALUE: optim MINIMIZES
	return(logLik)
}

