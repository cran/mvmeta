mvmeta.em <-
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

	#######################
	# STEP 1: EXPECTATION
	#######################
	
	# COMPUTE E(b) AND V(b) GIVEN Psi AND coef)
	Eblist <- mapply(function(invU,y,kX) {
		Psi%*%(tcrossprod(invU)%*%(y-kX%*%beta))},
		ls$invUlist,ylist,kXlist,SIMPLIFY=FALSE)
	Vblist <- lapply(ls$invUlist,function(invU) {
		Psi - Psi%*%(tcrossprod(invU)%*%Psi)})

	#######################
	# STEP 2: MAXIMIZATION
	#######################

	# COMPUTE E(b%*%t(b)
	Ebtblist <- mapply(function(Vb,Eb) Vb + tcrossprod(Eb),
		Vblist,Eblist,SIMPLIFY=FALSE)
        
	# NEW ESTIMATE FOR Psi
	Psi <- sumlist(Ebtblist)/m

	return(Psi)
}

