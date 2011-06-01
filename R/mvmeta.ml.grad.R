mvmeta.ml.grad <-
function(par, ylist, Slist, kXlist, nalist, nobs, k) {

	# RETRIEVE THE UPPER TRIANGULAR CHOLESKY MATRIX AND Psi
	L <- matrix(0,k,k)
	L[lower.tri(L,diag=TRUE)] <- par
	U <- t(L)
	Psi <- crossprod(U)

	# GET THE ESTIMATE OF beta CONDITIONAL ON Psi_theta
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

	# COMPUTE THE MATRIX DERIVATIVES OF EACH PARAMETER
	Dlist <- lapply(seq(length(par)),gradchol,k,U)

	# FUNCTION TO COMPUTE THE GRADIENT COMPONENTS FOR A GIVEN PARAMETER
	fgrad <- function(D,y,kX,invU,na) {
		sum(mapply(function(y,kX,invU,na) {
			invSigma <- tcrossprod(invU)
			res <- as.numeric(y-kX%*%beta)
			A <- crossprod(res,invSigma)%*%D[na,na]%*%invSigma%*%res
			B <- sum(diag(invSigma%*%D[na,na]))
			return(as.numeric(0.5*(A-B)))},y,kX,invU,na))
	}

	# REPEAT FOR ALL THE PARAMETERS 
	grad <- sapply(Dlist,fgrad,y=ylist,kX=kXlist,invU=invUlist,nalist)
	return(grad)
}

