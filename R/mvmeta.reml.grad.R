mvmeta.reml.grad <-
function(par, ylist, Slist, kXlist, nalist, nobs, k) {

	# RETRIEVE THE UPPER TRIANGULAR CHOLESKY MATRIX AND Psi
	L <- matrix(0,k,k)
	L[lower.tri(L,diag=TRUE)] <- par
	U <- t(L)
	Psi <- crossprod(U)

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

	tXMXtot <- sumlist(lapply(invtUXlist,function(x)crossprod(x)))
	invtXMXtot <- chol2inv(chol(tXMXtot))
	

	# COMPUTE THE MATRIX DERIVATIVES OF EACH PARAMETER
	Dlist <- lapply(seq(length(par)),gradchol,k,U)

	# FUNCTION TO COMPUTE THE GRADIENT COMPONENTS FOR A GIVEN PARAMETERS
	fgrad <- function(D,y,kX,invU,na) {
		sum(mapply(function(y,kX,invU,na) {
			invSigma <- tcrossprod(invU)
			MDM <- invSigma%*%D[na,na]%*%invSigma
			res <- as.numeric(y-kX%*%coef)
			A <- crossprod(res,MDM)%*%res
			B <- sum(diag(invSigma%*%D[na,na]))
			C <- sum(diag(invtXMXtot%*%crossprod(kX,MDM)%*%kX))
			return(as.numeric(0.5*(A-B+C)))},y,kX,invU,na))
	}

	# REPEAT FOR ALL THE PARAMETERS 
	grad <- sapply(Dlist,fgrad,y=ylist,kX=kXlist,invU=invUlist,nalist)
	return(grad)
}

